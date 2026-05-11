import { h } from "hastscript";
import addClasses from "rehype-class-names";
import rehypeAutolinkHeadings from "rehype-autolink-headings";
import rehypeKatex from "rehype-katex";
import rehypePrettyCode from "rehype-pretty-code";
import rehypeSlug from "rehype-slug";
import urls from "rehype-urls";
import rehypeWrap from "rehype-wrap-all";
import remarkDirective from "remark-directive";
import emoji from "remark-emoji";
import remarkGfm from "remark-gfm";
import remarkMath from "remark-math";
import type { AstroMarkdownOptions } from "@astrojs/markdown-remark";
import { transformerNotationDiff } from "@shikijs/transformers";
import admonitionsPlugin from "./remark-admonitions.js";
import { rehypeCheckboxParser } from "./rehype-checkbox-parser.ts";
import { rehypeHeadingNumbers } from "./rehype-heading-numbers.ts";
import { mermaid } from "./remark-mermaid.ts";

function urlTransformer(url: any) {
    const regex = /^https:\/\/(raw\.)*github/;
    if (!regex.test(url.href) && url.href?.endsWith(".md")) {
        url.href = url.href.replace(/\.md$/, "/");
        url.pathname = url.pathname.replace(/\.md$/, "/");
        url.path = url.path?.replace(/\.md$/, "/");
    } else if (!regex.test(url.href) && url.href?.endsWith(".mdx")) {
        url.href = url.href.replace(/\.mdx$/, "/");
        url.pathname = url.pathname.replace(/\.mdx$/, "/");
        url.path = url.path?.replace(/\.mdx$/, "/");
    }
}

export const sharedMarkdownConfig = {
    syntaxHighlight: false as const,
    remarkPlugins: [emoji, remarkGfm, remarkDirective, admonitionsPlugin, mermaid, remarkMath],
    rehypePlugins: [
        rehypeSlug,
        [rehypeAutolinkHeadings, { behavior: "append", content: h("i.ms-1.fas.fa-link.fa-xs.invisible") }],
        [addClasses, { table: "table table-hover table-sm small" }],
        [rehypeWrap, { selector: "table", wrapper: "div.table-responsive" }],
        [urls, urlTransformer],
        rehypeCheckboxParser,
        rehypeHeadingNumbers,
        [
            rehypePrettyCode,
            {
                defaultLang: "plaintext",
                keepBackground: true,
                theme: { dark: "github-dark-dimmed", light: "github-light" },
                transformers: [transformerNotationDiff()],
            },
        ],
        [rehypeKatex, { strict: false }],
    ],
} satisfies AstroMarkdownOptions;
