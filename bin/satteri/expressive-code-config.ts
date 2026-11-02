// Expressive Code (https://expressive-code.com/) config for satteri-rendered markdown.
//
// Shared by the `astro-expressive-code` integration in astro.config.base.mjs
// (Astro content), satteri-expressive-code in markdownConfig.ts (`renderMarkdown()`
// and `createSatteriPluginSets()` consumers like pipelineMarkdown.ts, which the
// integration can't reach), and hast-inline-code.ts (themes and language aliases).
//
// satteri-expressive-code still depends on expressive-code 0.42; the root
// package.json overrides it to our expressive-code version, so both paths run the
// same engine version with the plugins below.
import { pluginLineNumbers } from "@expressive-code/plugin-line-numbers";
import nextflowSvg from "../../sites/main-site/src/icons/logos/nextflow.svg?raw";
import { pluginCaptions } from "./ec-plugin-captions.ts";
import { pluginFileIcons } from "./ec-plugin-file-icons.ts";

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
 * astro-expressive-code integration and satteri-expressive-code.
 * A factory (not a constant) so each engine gets its own plugin instances.
 */
export function createEcConfig() {
    return {
        themes: [THEMES.dark, THEMES.light] as any,
        // The site toggles Bootstrap's data-bs-theme on <html> (ThemeSwitch.svelte);
        // the prefers-color-scheme media query covers pages before the attribute is set.
        themeCssSelector: (theme: { type: string }) => `[data-bs-theme='${theme.type}']`,
        // Line numbers stay opt-in via `showLineNumbers`, as with rehype-pretty-code.
        defaultProps: {
            showLineNumbers: false,
            // EC auto-detects a terminal frame for shell languages; opt out so
            // bash/sh blocks render as plain code frames. `overridesByLang` is a
            // `defaultProps` key, not a top-level EC option.
            overridesByLang: {
                "bash,sh": { frame: "none" },
            },
        },
        plugins: [
            pluginLineNumbers(),
            pluginCaptions(),
            pluginFileIcons([
                {
                    name: "nextflow",
                    svg: nextflowSvg,
                    extensions: [".nf", ".nf.test", ".config"],
                },
            ]),
        ],
        shiki: {
            langAlias: LANG_ALIAS,
        },
    };
}
