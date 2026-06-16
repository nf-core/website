// Expressive Code (https://expressive-code.com/) for satteri-rendered markdown.
//
// Astro content is handled by the `astro-expressive-code` integration in
// astro.config.base.mjs, which pushes its own satteri plugin into the shared
// processor (supported since astro-expressive-code 0.43.0). The plugin below
// covers the render paths the integration can't reach — `renderMarkdown()` and
// `createSatteriPluginSets()` consumers like pipelineMarkdown.ts — using the
// same engine and config, but inlining the engine's base/theme styles and JS
// modules once per document (like plain rehype-expressive-code usage).
//
// Use the exported factory directly in `hastPlugins` so the per-document
// style/script bookkeeping resets on every compile.
import { fileURLToPath } from "node:url";
import { pluginLineNumbers } from "@expressive-code/plugin-line-numbers";
import { createRenderer } from "rehype-expressive-code";
import { defineHastPlugin } from "satteri";
import { pluginCaptions } from "./ec-plugin-captions.ts";
import { pluginFileIcons } from "./ec-plugin-file-icons.ts";

const TAB_WIDTH = 2;

/** Shiki themes used for all code rendering, also by hast-inline-code.ts. */
export const THEMES = { dark: "github-dark-dimmed", light: "github-light" } as const;

/** Languages used in content that Shiki doesn't ship, also used by hast-inline-code.ts. */
export const LANG_ALIAS: Record<string, string> = {
    nextflow: "groovy",
    default: "text",
    none: "text",
    tree: "text",
};

/**
 * One source of truth for the Expressive Code config, shared between the
 * astro-expressive-code integration and the standalone satteri plugin below.
 * A factory (not a constant) so each engine gets its own plugin instances.
 */
export function createEcConfig() {
    return {
        themes: [THEMES.dark, THEMES.light] as any,
        // The site toggles Bootstrap's data-bs-theme on <html> (ThemeSwitch.svelte);
        // the prefers-color-scheme media query covers pages before the attribute is set.
        themeCssSelector: (theme: { type: string }) => `[data-bs-theme='${theme.type}']`,
        // Line numbers stay opt-in via `showLineNumbers`, as with rehype-pretty-code.
        defaultProps: { showLineNumbers: false },
        plugins: [
            pluginLineNumbers(),
            pluginCaptions(),
            pluginFileIcons([
                {
                    name: "nextflow",
                    path: fileURLToPath(
                        new URL("../../sites/main-site/src/icons/logos/nextflow.svg", import.meta.url),
                    ),
                    extensions: [".nf", ".nf.test", ".config"],
                },
            ]),
        ],
        shiki: {
            langAlias: LANG_ALIAS,
        },
    };
}

// One engine for all documents: theme loading is expensive and the output only
// depends on the block, not the document.
let rendererPromise: ReturnType<typeof createRenderer> | undefined;

function getRenderer() {
    rendererPromise ??= createRenderer(createEcConfig());
    return rendererPromise;
}

export function expressiveCodePlugin() {
    const addedStyles = new Set<string>();
    const addedJsModules = new Set<string>();

    return defineHastPlugin({
        name: "expressive-code",
        element: {
            filter: ["pre"],
            async visit(node: any, ctx) {
                // Only plain `pre > code > text` blocks; Expressive Code's own
                // rendered output (spans per line) never matches this shape.
                if (node.children?.length !== 1) return;
                const code = node.children[0];
                if (code?.type !== "element" || code.tagName !== "code") return;
                if (code.children?.length !== 1 || code.children[0]?.type !== "text") return;

                // Satteri keeps the fenced info-string language on data.lang
                // and the rest on data.meta.
                const lang: string = code.data?.lang ?? "";
                const meta: string = code.data?.meta ?? "";
                const value = ctx
                    .textContent(code)
                    .replace(/\n$/, "")
                    .replace(/\t/g, " ".repeat(TAB_WIDTH));

                const { ec, baseStyles, themeStyles, jsModules } = await getRenderer();
                const { renderedGroupAst, styles } = await ec.render({
                    code: value,
                    language: lang,
                    meta,
                });

                // Inline assets exactly once per document, like rehype-expressive-code.
                const extraElements: any[] = [];
                const stylesToPrepend: string[] = [];
                for (const style of [baseStyles, themeStyles, ...styles]) {
                    if (!style || addedStyles.has(style)) continue;
                    addedStyles.add(style);
                    stylesToPrepend.push(style);
                }
                if (stylesToPrepend.length) {
                    extraElements.push({
                        type: "element",
                        tagName: "style",
                        properties: {},
                        children: [{ type: "text", value: stylesToPrepend.join("") }],
                    });
                }
                for (const moduleCode of jsModules) {
                    if (addedJsModules.has(moduleCode)) continue;
                    addedJsModules.add(moduleCode);
                    extraElements.push({
                        type: "element",
                        tagName: "script",
                        properties: { type: "module" },
                        children: [{ type: "text", value: moduleCode }],
                    });
                }
                if (extraElements.length) {
                    const firstChild = renderedGroupAst.children[0] as any;
                    const firstChildIsStyle =
                        firstChild?.type === "element" && ["style", "link"].includes(firstChild?.tagName);
                    renderedGroupAst.children.splice(firstChildIsStyle ? 1 : 0, 0, ...extraElements);
                }
                return renderedGroupAst as any;
            },
        },
    });
}

export default expressiveCodePlugin;
