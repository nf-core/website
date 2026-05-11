import mdx from "@astrojs/mdx";
import netlify from "@astrojs/netlify";
import partytown from "@astrojs/partytown";
import sitemap from "@astrojs/sitemap";
import svelte from "@astrojs/svelte";
import yaml from "@rollup/plugin-yaml";
import { envField, fontProviders, svgoOptimizer } from "astro/config";
import markdownIntegration from "@mashehu/astropub-md";
import icon from "astro-icon";
import { sharedMarkdownConfig } from "./bin/markdownConfig.ts";

/**
 * Base Astro configuration shared across all nf-core subsites.
 * Import this and spread it into defineConfig() along with site-specific overrides.
 */
export default {
    site: "https://nf-co.re/",
    output: "static",
    adapter: netlify(),
    prefetch: false,
    env: {
        schema: {
            GITHUB_TOKEN: envField.string({
                context: "server",
                access: "secret",
                optional: false,
            }),
        },
    },
    fonts: [
        {
            provider: fontProviders.fontsource(),
            name: "Inter",
            cssVariable: "--font-inter",
            fallbacks: ["sans-serif"],
            weights: ["300 700"],
        },
        {
            provider: fontProviders.fontsource(),
            name: "Maven Pro",
            cssVariable: "--font-maven-pro",
            fallbacks: ["sans-serif"],
            weights: ["300 700"],
        },
    ],
    experimental: {
        svgOptimizer: svgoOptimizer(),
    },
    integrations: [
        svelte(),
        icon({
            include: {
                "file-icons": ["nextflow"],
                logos: [
                    "twitter",
                    "mastodon-icon",
                    "slack-icon",
                    "aws",
                    "microsoft-azure",
                    "github-actions",
                    "youtube-icon",
                    "linkedin",
                ],
                fa: ["github"],
                "fa-brands": ["github"],
                mdi: ["aws", "slack", "youtube", "hammer-wrench", "sprout"],
                octicon: [
                    "chevron-right-16",
                    "git-pull-request-16",
                    "law-16",
                    "link-external-16",
                    "mortar-board-16",
                    "play-16",
                    "table-16",
                    "tasklist-16",
                    "terminal-16",
                    "tools-16",
                    "question-16",
                    "book-16",
                ],
                "simple-icons": ["bluesky"],
            },
        }),
        sitemap(),
        partytown({
            config: {
                forward: ["dataLayer.push"],
            },
        }),
        mdx(),
        markdownIntegration(),
    ],
    build: {
        inlineStylesheets: "auto",
        format: "file",
    },
    vite: {
        plugins: [yaml()],
        ssr: {
            noExternal: ["@popperjs/core", "svelte-exmarkdown", "svelte-confetti", "@mashehu/astropub-md"],
        },
        resolve: {
            preserveSymlinks: true,
            browser: true,
            noExternal: ["@popperjs/core", "svelte-exmarkdown", "svelte-confetti", "@mashehu/astropub-md"],
        },
        css: {
            preprocessorOptions: {
                scss: {
                    silenceDeprecations: [
                        "legacy-js-api",
                        "color-functions",
                        "global-builtin",
                        "import",
                        "if-function",
                    ],
                },
            },
        },
    },
    image: {
        domains: [
            "raw.githubusercontent.com",
            "unsplash.com",
            "avatars.githubusercontent.com",
            "github.com",
            "nf-core-docs.netlify.app",
            "nf-core-main-site.netlify.app",
        ],
        service: {
            entrypoint: "astro/assets/services/sharp",
        },
    },
    markdown: sharedMarkdownConfig,
};
