// Satteri replacement for `rehype-katex` (paired with the built-in `math` feature,
// which replaces `remark-math`).
// Satteri parses $...$ / $$...$$ into mdast `inlineMath` / `math` nodes; this plugin
// renders them straight to HTML with KaTeX and splices the result in as raw HTML.
import katex from "katex";
import { defineMdastPlugin } from "satteri";

const render = (value: string, displayMode: boolean) =>
    katex.renderToString(value, { displayMode, strict: false, throwOnError: false });

export const mathKatexPlugin = defineMdastPlugin({
    name: "math-katex",
    math(node) {
        return { rawHtml: render(node.value, true) };
    },
    inlineMath(node) {
        // `{ rawHtml }` is spliced as flow content, which would wrap inline math in
        // its own <p>; an inline mdast html node keeps it inside the paragraph.
        return { type: "html", value: render(node.value, false) };
    },
});
