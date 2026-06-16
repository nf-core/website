import { createSatteriMarkdownProcessor } from "@astrojs/markdown-satteri";
import type { MarkdownHeading } from "astro";
import { createSatteriPluginSets } from "@root/bin/satteri/markdownConfig.ts";
import { createGitHubMarkdownPlugin, trimGitHubReadme } from "@root/bin/satteri/mdast-github-markdown.ts";

export async function renderPipelineMarkdown(
    body: string,
    repo: string,
    ref: string,
    parentDir: string,
): Promise<{ html: string; headings: MarkdownHeading[] }> {
    const { features, mdastPlugins, hastPlugins } = createSatteriPluginSets();
    // The processor is created per call: the GitHub URL rewriting is parameterised by
    // repo/ref, and the shared hast plugins carry per-document state. Creation is
    // cheap — satteri loads once and syntaxHighlight is off (hast-expressive-code.ts
    // renders the code blocks, with a shared Expressive Code engine).
    const processor = await createSatteriMarkdownProcessor({
        syntaxHighlight: false,
        features,
        mdastPlugins: [createGitHubMarkdownPlugin({ repo, ref, parentDirectory: parentDir }), ...mdastPlugins],
        hastPlugins,
    });
    const { code: html, metadata } = await processor.render(trimGitHubReadme(body));
    return { html, headings: (metadata?.headings ?? []) as MarkdownHeading[] };
}
