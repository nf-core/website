// Expressive Code plugin: support a `caption="..."` meta option that renders a
// caption below the code block, replicating rehype-pretty-code's captions
// (https://rehype-pretty.pages.dev/#captions).
//
// Usage in markdown:
//   ```groovy title="main.nf" caption="Versions emission block"
//   ...
//   ```
//
// EC parses `caption="..."` into `codeBlock.metaOptions` for free; this plugin
// reads it in `postprocessRenderedBlock` and appends a `<figcaption>` to the
// rendered `<figure>`. Registered in createEcConfig() so it applies to both the
// astro-expressive-code integration and the standalone satteri path.
import { definePlugin } from "@expressive-code/core";
import { h } from "hastscript";

export function pluginCaptions() {
    return definePlugin({
        name: "captions",
        hooks: {
            postprocessRenderedBlock: ({ codeBlock, renderData }) => {
                const caption = codeBlock.metaOptions.getString("caption");
                if (!caption) return;

                // renderData.blockAst is the block's root `<figure>` element.
                const figure = renderData.blockAst as any;
                figure.children ??= [];
                figure.children.push(h("figcaption", { className: ["ec-caption"] }, caption));
            },
        },
    });
}

export default pluginCaptions;
