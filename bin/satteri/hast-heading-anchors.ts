// Satteri replacement for `rehype-slug` + `rehype-autolink-headings`
// (behavior: "append", with the invisible link icon).
// Factory so the slugger's duplicate counter resets between documents.
import GithubSlugger from "github-slugger";
import { h } from "hastscript";
import { defineHastPlugin } from "satteri";

export function headingAnchorsPlugin() {
    const slugger = new GithubSlugger();
    return defineHastPlugin({
        name: "heading-anchors",
        element: {
            filter: ["h1", "h2", "h3", "h4", "h5", "h6"],
            visit(node: any, ctx) {
                let id = node.properties?.id;
                if (!id) {
                    id = slugger.slug(ctx.textContent(node));
                    ctx.setProperty(node, "id", id);
                }
                ctx.appendChild(
                    node,
                    h("a", { href: `#${id}` }, [h("i", { class: "ms-1 fas fa-link fa-xs invisible" })]),
                );
            },
        },
    });
}
