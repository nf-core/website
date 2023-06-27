import admonitionsPlugin from './bin/remark-admonitions.js';
import { mermaid } from './bin/remark-mermaid.ts';
import githubDarkDimmed from '/public/themes/github-dark-dimmed.json';
import githubLightTheme from '/public/themes/github-light.json';
import mdx from '@astrojs/mdx';
import netlify from '@astrojs/netlify/functions';
import partytown from '@astrojs/partytown';
import prefetch from '@astrojs/prefetch';
import sitemap from '@astrojs/sitemap';
import svelte from '@astrojs/svelte';
import markdownIntegration from '@astropub/md';
import yaml from '@rollup/plugin-yaml';
import { defineConfig } from 'astro/config';
import { FontaineTransform } from 'fontaine';
import { h } from 'hastscript';
import addClasses from 'rehype-add-classes';
import rehypeAutolinkHeadings from 'rehype-autolink-headings';
import rehypePrettyCode from 'rehype-pretty-code';
import rehypeSlug from 'rehype-slug';
import urls from 'rehype-urls';
import rehypeWrap from 'rehype-wrap-all';
import remarkDirective from 'remark-directive';
import emoji from 'remark-emoji';
import remarkGfm from 'remark-gfm';


// https://astro.build/config
export default defineConfig({
    site: 'https://astro--nf-core.netlify.app/',
    // TODO: switch back to 'https://nf-co.re/'
    output: 'hybrid',
    experimental: {
        assets: true,
    },
    adapter: netlify(),
    trailingSlash: 'always',
    integrations: [
        svelte(),
        sitemap(),
        markdownIntegration(),
        prefetch(),
        partytown({
            // Adds dataLayer.push as a forwarding-event.
            config: {
                forward: ['dataLayer.push'],
            },
        }),
        mdx(),
    ],
    build: {
        inlineStylesheets: 'auto',
    },
    vite: {
        plugins: [
            yaml(),
            FontaineTransform.vite({ // avoid flash of unstyled text by interjecting fallback system fonts
                fallbacks: ['BlinkMacSystemFont', 'Segoe UI', 'Helvetica Neue', 'Arial', 'Noto Sans'],
            }),
        ],
        ssr: {
            noExternal: ['@popperjs/core', 'bin/cache.js'],
        },
        resolve: {
            preserveSymlinks: true,
        },
    },
    image: {
        service: {
            entrypoint: 'astro/assets/services/sharp',
        },
    },
    markdown: {
        syntaxHighlight: false,
        remarkPlugins: [emoji, remarkGfm, remarkDirective, admonitionsPlugin, mermaid],
        // NOTE: Also update the plugins in `src/components/Markdown.svelte`!
        rehypePlugins: [
            rehypeSlug,
            [
                rehypeAutolinkHeadings,
                {
                    behavior: 'append',
                    content: h('i.ms-1.fas.fa-link.invisible'),
                },
            ],
            [
                addClasses,
                {
                    table: 'table table-hover table-sm small',
                },
            ],
            [
                rehypeWrap,
                {
                    selector: 'table',
                    wrapper: 'div.table-responsive',
                },
            ],
            [
                urls,
                (url) => {
                    if (url.href?.endsWith('.md')) {
                        url.href = url.href.replace(/\.md$/, '/');
                        url.pathname = url.pathname.replace(/\.md$/, '/');
                        url.path = url.path.replace(/\.md$/, '/');
                    }
                },
            ],
            [
                rehypePrettyCode,
                {
                    langPrefix: 'language-',
                    keepBackground: true,
                    theme: {
                        dark: 'github-dark-dimmed',
                        light: githubLightTheme,
                    },
                },
            ],
        ],
    },
});