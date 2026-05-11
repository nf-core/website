import { createMarkdownProcessor } from "@astrojs/markdown-remark";
import type { MarkdownHeading } from "astro";
import { sharedMarkdownConfig } from "@root/bin/markdownConfig.ts";
import remarkGitHubMarkdown from "@root/bin/remark-github-markdown.js";

// Created once for the entire build — Shiki grammars and plugin instances are reused.
const processorPromise = createMarkdownProcessor({
    ...sharedMarkdownConfig,
    remarkPlugins: [remarkGitHubMarkdown, ...sharedMarkdownConfig.remarkPlugins],
});

export async function renderPipelineMarkdown(
    body: string,
    repo: string,
    ref: string,
    parentDir: string,
): Promise<{ html: string; headings: MarkdownHeading[] }> {
    const processor = await processorPromise;
    const { code: html, metadata } = await processor.render(body, {
        frontmatter: { repo, ref, parent_directory: parentDir },
    });
    return { html, headings: (metadata?.headings ?? []) as MarkdownHeading[] };
}
