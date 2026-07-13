// Expressive Code plugin: inline a custom SVG logo as the frame icon for code
// blocks whose title matches given file-name suffixes (e.g. the Nextflow logo
// for `.nf`/`.config` files).
//
// Modelled on https://github.com/xt0rted/expressive-code-file-icons: it hooks
// `postprocessRenderedBlock`, finds the frame's `<span class="title">` and
// prepends the inline logo. Registered once in createEcConfig() so it applies
// to both the astro-expressive-code integration and the standalone satteri path.
import { definePlugin } from "@expressive-code/core";
import { fromHtml } from "hast-util-from-html";

export interface FileIcon {
    /** Identifier used for the injected SVG's CSS class (`<name>-file-icon`). */
    name: string;
    /** SVG file contents, inlined at bundle time (e.g. via `import svg from "...svg?raw"`). */
    svg: string;
    /** File-name suffixes that map to this icon, e.g. `[".nf", ".config"]`. */
    extensions: string[];
}

type LoadedIcon = {
    extensions: string[];
    className: string;
    template: any;
};

function loadIcon({ name, svg, extensions }: FileIcon): LoadedIcon {
    const template = fromHtml(svg, { fragment: true, space: "svg" }).children.find(
        (node: any) => node.type === "element",
    );
    return { extensions, className: `${name}-file-icon`, template };
}

// Recursively find the frame's title span (`<span class="title">`).
function findTitleSpan(node: any): any {
    if (!node || typeof node !== "object") return undefined;
    if (node.type === "element" && node.tagName === "span") {
        const classes = node.properties?.className;
        const classList = Array.isArray(classes) ? classes : typeof classes === "string" ? classes.split(/\s+/) : [];
        if (classList.includes("title")) return node;
    }
    for (const child of node.children ?? []) {
        const found = findTitleSpan(child);
        if (found) return found;
    }
    return undefined;
}

export function pluginFileIcons(icons: FileIcon[]) {
    // Parse each SVG once; every insertion gets a clone (hast nodes are mutated
    // downstream, so the template must never be shared).
    const loaded = icons.map(loadIcon);

    return definePlugin({
        name: "file-icons",
        hooks: {
            postprocessRenderedBlock: ({ codeBlock, renderData }) => {
                const title = codeBlock.props.title;
                if (!title) return;

                const match = loaded.find((icon) => icon.extensions.some((ext) => title.endsWith(ext)));
                if (!match?.template) return;

                const titleSpan = findTitleSpan(renderData.blockAst);
                if (!titleSpan) return;

                const node = structuredClone(match.template);
                node.properties ??= {};
                const className = node.properties.className;
                const existing = Array.isArray(className) ? className : [];
                // Bootstrap utilities for spacing/alignment, matching the runtime
                // file-icons in MarkdownWrapper.svelte.
                node.properties.className = [...existing, match.className, "me-2", "align-middle"];
                node.properties["aria-hidden"] = "true";

                titleSpan.children ??= [];
                titleSpan.children.unshift(node);
            },
        },
    });
}

export default pluginFileIcons;
