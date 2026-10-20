// Satteri replacement for `rehype-urls` with the urlTransformer from
// bin/markdownConfig.ts: rewrite non-GitHub links to .md/.mdx files to
// directory-style URLs ("page.md" → "page/", "page.md#section" → "page/#section").
import { defineHastPlugin } from "satteri";

const GITHUB_URL = /^https:\/\/(raw\.)*github/;
// the .md/.mdx extension, followed by the end of the link, an #anchor or a ?query
const MARKDOWN_FILE = /\.mdx?(?=[#?]|$)/;
// a path ending in .md/.mdx, followed by the end of the link, an #anchor or a ?query
const MARKDOWN_FILE = /^([^?#]*)\.mdx?(?=[?#]|$)/;

export const markdownLinksPlugin = defineHastPlugin({
    name: "markdown-links",
    element: {
        filter: ["a"],
        visit(node: any, ctx) {
            const href = node.properties?.href;
            if (typeof href !== "string" || GITHUB_URL.test(href)) return;
            if (MARKDOWN_FILE.test(href)) {
                ctx.setProperty(node, "href", href.replace(MARKDOWN_FILE, "/"));
            }
        },
    },
});
