import admonitionsPlugin from "./bin/remark-admonitions";
import { mermaid } from "./bin/remark-mermaid.ts";
import { rehypeCheckboxParser } from "./bin/rehype-checkbox-parser.ts";
import { rehypeHeadingNumbers } from "./bin/rehype-heading-numbers.ts";
import mdx from "@astrojs/mdx";
import netlify from "@astrojs/netlify";
import partytown from "@astrojs/partytown";
import sitemap from "@astrojs/sitemap";
import svelte from "@astrojs/svelte";
import yaml from "@rollup/plugin-yaml";
import { envField, fontProviders } from "astro/config";
import { h } from "hastscript";
import addClasses from "rehype-class-names";
import rehypeAutolinkHeadings from "rehype-autolink-headings";
import rehypeKatex from "rehype-katex";
import rehypePrettyCode from "rehype-pretty-code";
import { transformerNotationDiff } from "@shikijs/transformers";
import rehypeSlug from "rehype-slug";
import urls from "rehype-urls";
import rehypeWrap from "rehype-wrap-all";
import remarkDirective from "remark-directive";
import emoji from "remark-emoji";
import remarkGfm from "remark-gfm";
import remarkMath from "remark-math";
import markdownIntegration from "@mashehu/astropub-md";
import icon from "astro-icon";

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
    experimental: {
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
        svgo: true,
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
                mdi: ["aws", "slack", "youtube"],
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
            noExternal: ["@popperjs/core"],
        },
        resolve: {
            preserveSymlinks: true,
            browser: true,
        },
        css: {
            preprocessorOptions: {
                scss: {
                    api: "modern-compiler",
                    silenceDeprecations: ["legacy-js-api", "color-functions", "global-builtin", "import"],
                },
            },
        },
    },
    image: {
        domains: ["raw.githubusercontent.com", "unsplash.com", "avatars.githubusercontent.com", "github.com"],
        service: {
            entrypoint: "astro/assets/services/sharp",
        },
    },
    markdown: {
        syntaxHighlight: false,
        remarkPlugins: [emoji, remarkGfm, remarkDirective, admonitionsPlugin, mermaid, remarkMath],
        rehypePlugins: [
            rehypeSlug,
            [
                rehypeAutolinkHeadings,
                {
                    behavior: "append",
                    content: h("i.ms-1.fas.fa-link.fa-xs.invisible"),
                },
            ],
            [
                addClasses,
                {
                    table: "table table-hover table-sm small",
                },
            ],
            [
                rehypeWrap,
                {
                    selector: "table",
                    wrapper: "div.table-responsive",
                },
            ],
            [
                urls,
                (url) => {
                    const regex = /^https:\/\/(raw.)*github/;
                    if (!regex.test(url.href) && url.href?.endsWith(".md")) {
                        url.href = url.href.replace(/\.md$/, "/");
                        url.pathname = url.pathname.replace(/\.md$/, "/");
                        url.path = url.path.replace(/\.md$/, "/");
                    } else if (!regex.test(url.href) && url.href?.endsWith(".mdx")) {
                        url.href = url.href.replace(/\.mdx$/, "/");
                        url.pathname = url.pathname.replace(/\.mdx$/, "/");
                        url.path = url.path.replace(/\.mdx$/, "/");
                    }
                },
            ],
            rehypeCheckboxParser,
            rehypeHeadingNumbers,
            [
                rehypePrettyCode,
                {
                    defaultLang: "plaintext",
                    keepBackground: true,
                    theme: {
                        dark: "github-dark-dimmed",
                        light: "github-light",
                    },
                    transformers: [transformerNotationDiff()],
                },
            ],
            rehypeKatex,
        ],
    },
};
