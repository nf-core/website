// Inline code highlighting, restoring rehype-pretty-code's `code{:lang}` syntax
// (https://rehype-pretty.pages.dev/#inline-code), which is widely used in content.
// Expressive Code only renders fenced blocks, so inline code is highlighted with
// Shiki directly, using the same themes. Tokens carry --shiki-light/--shiki-dark
// variables (no inline color); main.scss switches the variables per theme.
import { defineHastPlugin } from "satteri";
import { codeToHast } from "shiki";
import { LANG_ALIAS, THEMES } from "./hast-expressive-code.ts";

const INLINE_LANG = /\{:([\w.+#-]+)\}\s*$/;

export const inlineCodePlugin = defineHastPlugin({
    name: "inline-code-shiki",
    element: {
        filter: ["code"],
        async visit(node: any, ctx) {
            // Fenced code blocks carry their language on data.lang and belong to
            // Expressive Code; inline code is a code element with a single text
            // child and no newlines, ending in a {:lang} marker.
            if (node.data?.lang !== undefined) return;
            if (node.children?.length !== 1 || node.children[0]?.type !== "text") return;
            const text: string = node.children[0].value;
            if (text.includes("\n")) return;
            const match = text.match(INLINE_LANG);
            if (!match) return;
            const value = text.slice(0, match.index);
            const lang = LANG_ALIAS[match[1]] ?? match[1];

            let inline: any;
            try {
                inline = await codeToHast(value, {
                    lang,
                    themes: THEMES,
                    defaultColor: false,
                    cssVariablePrefix: "--shiki-",
                    structure: "inline",
                });
            } catch {
                // Unknown language: keep the code as plain text, but drop the marker.
                return {
                    type: "element",
                    tagName: "code",
                    properties: { ...node.properties },
                    children: [{ type: "text", value }],
                };
            }
            return {
                type: "element",
                tagName: "code",
                properties: { ...node.properties, "data-language": match[1], "data-inline-code": "" },
                children: inline.children,
            };
        },
    },
});

export default inlineCodePlugin;
