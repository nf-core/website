// Satteri replacement for `remark-emoji`.
// Converts :shortcode: emojis in text nodes to their unicode equivalent.
import { emojify } from "node-emoji";
import { defineMdastPlugin } from "satteri";

export const emojiPlugin = defineMdastPlugin({
    name: "emoji",
    text(node, ctx) {
        const value = emojify(node.value);
        if (value !== node.value) {
            ctx.setProperty(node, "value", value);
        }
    },
});

export default emojiPlugin;
