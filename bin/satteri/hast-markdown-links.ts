// Satteri replacement for `rehype-urls` with the urlTransformer from
// bin/markdownConfig.ts: strip the .md/.mdx extension from non-GitHub links to
// markdown files ("page.md" → "page", "page.md#section" → "page#section"), since
// pages are built with `build.format: "file"` and served without a trailing slash.
import { defineHastPlugin } from "satteri";

const GITHUB_URL = /^https:\/\/(raw\.)*github/;
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
                ctx.setProperty(node, "href", href.replace(MARKDOWN_FILE, "$1"));
            }
        },
    },
});
