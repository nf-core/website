// Satteri replacement for bin/remark-mermaid.ts.
// Replaces ```mermaid code blocks with a placeholder div that is hydrated client-side.
import { defineMdastPlugin } from "satteri";

const escapeMap: Record<string, string> = {
    "&": "&amp;",
    "<": "&lt;",
    ">": "&gt;",
    '"': "&quot;",
    "'": "&#39;",
    // newlines in the mermaid source must be escaped too: the raw HTML has to stay
    // a single line, otherwise the attribute is split across markdown HTML blocks
    "\n": "&#10;",
    "\r": "&#13;",
};

const escapeHtml = (str: string) => str.replace(/[&<>"'\n\r]/g, (c) => escapeMap[c]);

export const mermaidPlugin = defineMdastPlugin({
    name: "mermaid",
    code(node) {
        if (node.lang !== "mermaid") return;

        return {
            rawHtml: `<div class="mermaid" data-content="${escapeHtml(node.value)}"><i class="mt-5 m-auto text-success fa-regular fa-spinner-third fa-spin fa-3x"></i><p>Loading graph</p></div>`,
        };
    },
});

export default mermaidPlugin;
