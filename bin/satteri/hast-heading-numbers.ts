// Satteri replacement for bin/rehype-heading-numbers.ts.
// Prepends section numbers to headings when the document frontmatter opts in via
// markdownPlugin: addNumbersToHeadings.
//
// The old plugin had a "minDepth" scan pass, but it initialised minDepth to 1 and
// only ever took Math.min, so it always resolved to 1 — the scan is dropped here
// and depth is simply the heading level, which preserves the old output exactly.
import { defineHastPlugin } from "satteri";
import type { MarkdownFlags } from "./frontmatter-flags.ts";

export function createHeadingNumbersPlugin(flags: MarkdownFlags) {
    return () => {
        let counters: number[] = [1];
        return defineHastPlugin({
            name: "heading-numbers",
            element: {
                filter: ["h1", "h2", "h3", "h4", "h5", "h6"],
                visit(node: any) {
                    if (!flags.headingNumbers) return;

                    const depth = parseInt(node.tagName.slice(1));

                    // Initialize counters for new depths
                    if (depth > counters.length) {
                        counters.push(1);
                    } else {
                        // Increment counter for current depth and reset counters for bigger depths
                        counters = counters.slice(0, depth);
                        counters[depth - 1] = counters[depth - 1] + 1;
                    }

                    const counterText = counters.slice(1, depth).join(".");

                    const children = [...(node.children ?? [])];
                    // remove hard-coded heading numbers
                    if (children[0]?.type === "text") {
                        children[0] = {
                            ...children[0],
                            value: children[0].value.replace(/^\s*\d+(\.\d+)*\s/, ""),
                        };
                    }
                    children.unshift({ type: "text", value: `${counterText} ` });

                    return { ...node, children };
                },
            },
        });
    };
}
