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
            ctx.setProperty(node, "className", TABLE_CLASSES);
            ctx.wrapNode(node, h("div", { class: "table-responsive" }));
        },
    },
});
