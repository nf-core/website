// Satteri port of bin/markdownConfig.ts.
//
// Plugin mapping from the old unified pipeline:
//   remark-gfm, remark-math, remark-directive  → built-in `features` (no plugin needed)
//   remark-emoji                               → mdast-emoji.ts
//   bin/remark-admonitions.js                  → admonitions.ts (mdast + hast halves)
//   bin/remark-mermaid.ts                      → mdast-mermaid.ts
//   rehype-slug + rehype-autolink-headings     → hast-heading-anchors.ts
//   rehype-class-names + rehype-wrap-all       → hast-tables.ts
//   rehype-urls (urlTransformer)               → hast-markdown-links.ts
//   bin/rehype-checkbox-parser.ts              → hast-checkbox-parser.ts
//   bin/rehype-heading-numbers.ts              → hast-heading-numbers.ts
//   rehype-pretty-code                         → hast-pretty-code.ts
//   rehype-katex                               → mdast-math-katex.ts
//
// `syntaxHighlight: false` keeps @astrojs/markdown-satteri's own highlight plugin
// out of the pipeline (hast-pretty-code.ts does the highlighting), mirroring the old
// config. Astro's wrapper appends its image plugins and a heading-ids plugin after
// our hast plugins; the heading-ids plugin respects the ids hast-heading-anchors.ts
// already set, so the `headings` metadata (TOCs) stays in sync.
import { satteri } from "@astrojs/markdown-satteri";
import { markdownToHtml } from "satteri";
import { admonitionsHastPlugin, admonitionsMdastPlugin, directiveTextMergePlugin } from "./admonitions.ts";
import { createMarkdownFlags, type MarkdownFlags } from "./frontmatter-flags.ts";
import { createCheckboxParser, type CheckboxEntry } from "./hast-checkbox-parser.ts";
import { createHeadingNumbersPlugin } from "./hast-heading-numbers.ts";
import headingAnchorsPlugin from "./hast-heading-anchors.ts";
import markdownLinksPlugin from "./hast-markdown-links.ts";
import prettyCodePlugin from "./hast-pretty-code.ts";
import tablesPlugin from "./hast-tables.ts";
import emojiPlugin from "./mdast-emoji.ts";
import mathKatexPlugin from "./mdast-math-katex.ts";
import mermaidPlugin from "./mdast-mermaid.ts";

const features = {
    math: true,
    directive: true,
    // gfm and smartPunctuation are set by @astrojs/markdown-satteri from the Astro
    // markdown config defaults (both on), matching the old unified pipeline.
};

const sharedMdastPlugins = [
    emojiPlugin,
    directiveTextMergePlugin,
    admonitionsMdastPlugin,
    mermaidPlugin,
    mathKatexPlugin,
] as any[];

const hastPlugins = (flags: MarkdownFlags, checkboxPlugin: unknown) =>
    [
        admonitionsHastPlugin,
        headingAnchorsPlugin,
        tablesPlugin,
        markdownLinksPlugin,
        checkboxPlugin,
        createHeadingNumbersPlugin(flags),
        prettyCodePlugin,
    ] as any[];

/**
 * A fresh, self-contained set of the shared satteri plugins plus the per-document
 * state they share. Use this to build additional processors (e.g. pipelineMarkdown)
 * that should render identically to the site-wide config.
 */
export function createSatteriPluginSets() {
    const flags: MarkdownFlags = { checklist: false, headingNumbers: false };
    const checkboxParser = createCheckboxParser(flags);
    return {
        features,
        mdastPlugins: [...sharedMdastPlugins],
        hastPlugins: hastPlugins(flags, checkboxParser.plugin),
        flags,
        checkboxParser,
    };
}

/**
 * Astro `markdown` config (drop-in replacement for the old sharedMarkdownConfig).
 *
 * The checkbox parser and heading numbers are gated on frontmatter (markdownPlugin:
 * checklist / addNumbersToHeadings). Astro parses frontmatter before rendering and
 * passes it per render — satteri plugins never see it — so the renderer is wrapped
 * to copy the flags into the shared plugin state before each render and to write the
 * collected checkboxes back to the frontmatter afterwards (the old plugins did the
 * same through the vfile). Renders are serialised because that state is shared.
 */
export function satteriSharedMarkdownConfig() {
    const { flags, checkboxParser, mdastPlugins, hastPlugins } = createSatteriPluginSets();

    const processor = satteri({
        features,
        mdastPlugins,
        hastPlugins,
    });

    const createRenderer = processor.createRenderer.bind(processor);
    processor.createRenderer = async (shared: any) => {
        const renderer = await createRenderer(shared);
        let queue: Promise<unknown> = Promise.resolve();
        return {
            ...renderer,
            render(content: string, renderOpts: any) {
                const result = queue.then(async () => {
                    const frontmatter = renderOpts?.frontmatter ?? {};
                    flags.checklist = !!frontmatter.markdownPlugin?.includes("checklist");
                    flags.headingNumbers = !!frontmatter.markdownPlugin?.includes("addNumbersToHeadings");
                    const rendered = await renderer.render(content, renderOpts);
                    if (flags.checklist) {
                        frontmatter.checkboxes = checkboxParser.getCheckboxes();
                    }
                    return rendered;
                });
                queue = result.then(
                    () => {},
                    () => {},
                );
                return result;
            },
        };
    };

    return {
        syntaxHighlight: false as const,
        processor,
    };
}

export interface RenderedMarkdown {
    html: string;
    frontmatter: unknown;
    checkboxes: CheckboxEntry[];
}

/**
 * Standalone rendering outside Astro (frontmatter stays in the source, so the
 * markdownPlugin flags are read from the mdast yaml node instead). Documents must
 * be rendered sequentially — the per-document state is shared between the plugin
 * factories.
 */
export async function renderMarkdown(source: string, filename?: string): Promise<RenderedMarkdown> {
    const { flags, plugin: frontmatterFlagsPlugin } = createMarkdownFlags();
    const checkboxParser = createCheckboxParser(flags);

    const { html, frontmatter } = await markdownToHtml(source, {
        features: { ...features, gfm: true, frontmatter: true },
        mdastPlugins: [frontmatterFlagsPlugin, ...sharedMdastPlugins],
        hastPlugins: hastPlugins(flags, checkboxParser.plugin),
        filename,
    });
    return { html, frontmatter, checkboxes: checkboxParser.getCheckboxes() };
}
