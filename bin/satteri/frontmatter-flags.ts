// The old rehype plugins read per-document opt-ins from the Astro vfile
// (file.data.astro.frontmatter.markdownPlugin). Satteri has no vfile, but the
// `frontmatter` feature exposes the YAML block as an mdast `yaml` node, and mdast
// plugins always run before hast plugins — so this plugin extracts the flags into a
// shared object the hast plugins read.
//
// Plugin factories are re-invoked for every document, which is what resets the flags
// between documents. Note: this assumes documents are processed sequentially.
import { defineMdastPlugin } from "satteri";

export interface MarkdownFlags {
    checklist: boolean;
    headingNumbers: boolean;
}

export function createMarkdownFlags() {
    const flags: MarkdownFlags = { checklist: false, headingNumbers: false };

    const plugin = () => {
        flags.checklist = false;
        flags.headingNumbers = false;
        return defineMdastPlugin({
            name: "frontmatter-flags",
            yaml(node) {
                // Grab the value of the `markdownPlugin` key (scalar, flow or block list).
                const block = node.value.match(/^markdownPlugin:((?:.*)(?:\r?\n[ \t-].*)*)/m)?.[1] ?? "";
                flags.checklist = block.includes("checklist");
                flags.headingNumbers = block.includes("addNumbersToHeadings");
            },
        });
    };

    return { flags, plugin };
}
