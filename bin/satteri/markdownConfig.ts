// Satteri port of bin/markdownConfig.ts.
//
// Plugin mapping from the old unified pipeline:
//   remark-gfm, remark-math, remark-directive  → built-in `features` (no plugin needed)
//   remark-emoji                               → satteri-emoji
//   bin/remark-admonitions.js                  → admonitions.ts (mdast + hast halves)
//   bin/remark-mermaid.ts                      → mdast-mermaid.ts
//   rehype-slug + rehype-autolink-headings     → satteri-slug + satteri-autolink-headings
//   rehype-class-names + rehype-wrap-all       → hast-tables.ts
//   rehype-urls (urlTransformer)               → hast-markdown-links.ts
//   bin/rehype-checkbox-parser.ts              → hast-checkbox-parser.ts
//   bin/rehype-heading-numbers.ts              → hast-heading-numbers.ts
//   rehype-katex                               → satteri-katex
//
// Code blocks are rendered by Expressive Code: the astro-expressive-code
// integration (astro.config.base.mjs) pushes its satteri plugin into the shared
// processor's options, and the standalone paths below (createSatteriPluginSets,
// renderMarkdown) append satteri-expressive-code with the same config
// (expressive-code-config.ts).
// `syntaxHighlight: false` keeps @astrojs/markdown-satteri's own highlight plugin
// out of the pipeline. Astro's wrapper appends its image plugins and a heading-ids plugin after
// our hast plugins; the heading-ids plugin respects the ids satteri-slug
// already set, so the `headings` metadata (TOCs) stays in sync.
import { createSatteriMarkdownProcessor, satteri } from "@astrojs/markdown-satteri";
import { h } from "hastscript";
import { markdownToHtml } from "satteri";
import { satteriAutolinkHeadings } from "satteri-autolink-headings";
import { satteriEmoji } from "satteri-emoji";
import expressiveCode from "satteri-expressive-code";
import { satteriKatex } from "satteri-katex";
import { satteriSlug } from "satteri-slug";
import { admonitionsHastPlugin, admonitionsMdastPlugin, directiveTextMergePlugin } from "./admonitions.ts";
import { createEcConfig } from "./expressive-code-config.ts";
import { createCheckboxParser } from "./hast-checkbox-parser.ts";
import { headingNumbersPlugin } from "./hast-heading-numbers.ts";
import { inlineCodePlugin } from "./hast-inline-code.ts";
import { markdownLinksPlugin } from "./hast-markdown-links.ts";
import { tablesPlugin } from "./hast-tables.ts";
import { mermaidPlugin } from "./mdast-mermaid.ts";

const features = {
    math: true,
    directive: true,
    // gfm and smartPunctuation are set by @astrojs/markdown-satteri from the Astro
    // markdown config defaults (both on), matching the old unified pipeline.
};

// A plugin factory: satteri calls it per document, so the once-per-document
// style/script bookkeeping resets, while the engine (theme loading is expensive)
// is created lazily once and shared by all documents.
const expressiveCodePlugin = expressiveCode(createEcConfig());

const sharedMdastPlugins = [
    // GitHub shortcode names (:tada:, :heart:, …), as used by node-emoji/remark-emoji before
    satteriEmoji({ shortcodes: ["github"] }),
    directiveTextMergePlugin,
    admonitionsMdastPlugin,
    mermaidPlugin,
    satteriKatex(),
] as any[];

// The checkbox parser and heading numbers are opt-in per document (markdownPlugin:
// checklist / addNumbersToHeadings in the frontmatter); the opt-in renderer below
// appends them to the pipeline of documents that ask for them.
// inlineCodePlugin runs before satteriSlug so the `{:lang}` markers of inline code
// in headings are gone before the heading text is slugged.
const buildHastPlugins = () =>
    [
        admonitionsHastPlugin,
        inlineCodePlugin,
        // A factory, so the slugger's duplicate counter resets between documents.
        satteriSlug(),
        satteriAutolinkHeadings({
            behavior: "append",
            content: () => h("i", { class: "ms-1 fas fa-link fa-xs invisible" }),
        }),
        tablesPlugin,
        markdownLinksPlugin,
    ] as any[];

/**
 * A fresh set of the shared satteri plugins. Use this to build additional processors
 * (e.g. pipelineMarkdown) that should render identically to the site-wide config.
 */
export function createSatteriPluginSets() {
    return {
        features,
        mdastPlugins: [...sharedMdastPlugins],
        // The astro-expressive-code integration only reaches the shared Astro
        // processor, so these stateless plugin sets bring their own copy.
        hastPlugins: [...buildHastPlugins(), expressiveCodePlugin],
    };
}

/**
 * Astro `markdown` config (drop-in replacement for the old sharedMarkdownConfig).
 *
 * The checkbox parser and heading numbers are gated on frontmatter (markdownPlugin:
 * checklist / addNumbersToHeadings). Astro parses frontmatter before rendering and
 * passes it per render — satteri plugins never see it — so the renderer routes
 * opt-in documents through their own per-render processor that includes the gated
 * plugins and writes the collected checkboxes back to the frontmatter (the old
 * plugins did the same through the vfile). Everything else goes through the shared
 * stateless renderer, so renders stay concurrent.
 */
export function satteriSharedMarkdownConfig() {
    const processor = satteri({
        features,
        mdastPlugins: sharedMdastPlugins,
        hastPlugins: buildHastPlugins(),
    });

    const createRenderer = processor.createRenderer.bind(processor);
    processor.createRenderer = async (shared: any) => {
        const renderer = await createRenderer(shared);
        return {
            ...renderer,
            async render(content: string, renderOpts: any) {
                const frontmatter = renderOpts?.frontmatter ?? {};
                const checklist = !!frontmatter.markdownPlugin?.includes("checklist");
                const headingNumbers = !!frontmatter.markdownPlugin?.includes("addNumbersToHeadings");
                if (!checklist && !headingNumbers) {
                    return renderer.render(content, renderOpts);
                }

                const checkboxParser = checklist ? createCheckboxParser() : undefined;
                // Read the plugin arrays from processor.options (not the module
                // constants): integrations like astro-expressive-code push their
                // own plugins into them during astro:config:setup.
                const optInRenderer = await createSatteriMarkdownProcessor({
                    ...shared,
                    features: processor.options.features,
                    mdastPlugins: processor.options.mdastPlugins,
                    hastPlugins: [
                        ...processor.options.hastPlugins,
                        ...(checkboxParser ? [checkboxParser.plugin] : []),
                        ...(headingNumbers ? [headingNumbersPlugin] : []),
                    ],
                });
                const rendered = await optInRenderer.render(content, renderOpts);
                if (checkboxParser) {
                    frontmatter.checkboxes = checkboxParser.getCheckboxes();
                }
                return rendered;
            },
        };
    };

    return {
        syntaxHighlight: false as const,
        processor,
    };
}

/**
 * Standalone rendering outside Astro, for one-off markdown strings (e.g. pipeline
 * descriptions). gfm and smartPunctuation are filled in explicitly because raw
 * markdownToHtml defaults them off, unlike the Astro path.
 */
export async function renderMarkdown(source: string, filename?: string): Promise<{ html: string }> {
    const { html } = await markdownToHtml(source, {
        features: { ...features, gfm: true, smartPunctuation: true, frontmatter: true },
        mdastPlugins: sharedMdastPlugins,
        hastPlugins: [...buildHastPlugins(), expressiveCodePlugin],
        filename,
    });
    return { html };
}
