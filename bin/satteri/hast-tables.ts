// Satteri replacement for `rehype-class-names` ({ table: ... }) and
// `rehype-wrap-all` ({ selector: "table", wrapper: "div.table-responsive" }).
import { h } from "hastscript";
import { defineHastPlugin } from "satteri";

const TABLE_CLASSES = ["table", "table-hover", "table-sm", "small"];

export const tablesPlugin = defineHastPlugin({
    name: "tables",
    element: {
        filter: ["table"],
        visit(node: any, ctx) {
            const existing = node.properties?.className;
            const classes = Array.isArray(existing) ? existing : typeof existing === "string" ? existing.split(" ") : [];
            ctx.setProperty(node, "className", [...classes, ...TABLE_CLASSES.filter((c) => !classes.includes(c))]);
            ctx.wrapNode(node, h("div", { class: "table-responsive" }));
        },
    },
});

export default tablesPlugin;
