// Satteri replacement for `rehype-urls` with the urlTransformer from
// bin/markdownConfig.ts: rewrite non-GitHub links ending in .md/.mdx to
// directory-style URLs ("page.md" → "page/").
import { defineHastPlugin } from "satteri";

const GITHUB_URL = /^https:\/\/(raw\.)*github/;

export const markdownLinksPlugin = defineHastPlugin({
    name: "markdown-links",
    element: {
        filter: ["a"],
        visit(node: any, ctx) {
            const href = node.properties?.href;
            if (typeof href !== "string" || GITHUB_URL.test(href)) return;
            if (href.endsWith(".md") || href.endsWith(".mdx")) {
                ctx.setProperty(node, "href", href.replace(/\.mdx?$/, "/"));
            }
        },
    },
});

export default markdownLinksPlugin;
