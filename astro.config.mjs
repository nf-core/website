import admonitionsPlugin from './bin/remark-admonitions.js';
import { mermaid } from './bin/remark-mermaid.ts';
import pipelines_json from '/public/pipelines.json';
import githubDarkDimmed from '/public/themes/github-dark-dimmed.json';
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
import rehypeKatex from 'rehype-katex';
import rehypePrettyCode from 'rehype-pretty-code';
import rehypeSlug from 'rehype-slug';
import urls from 'rehype-urls';
import rehypeWrap from 'rehype-wrap-all';
import remarkDirective from 'remark-directive';
import emoji from 'remark-emoji';
import remarkGfm from 'remark-gfm';
import remarkMath from 'remark-math';


const latestToolsRelease = await fetch('https://api.github.com/repos/nf-core/tools/releases/latest')
    .then((res) => res.json())
    .then((json) => json.tag_name);
let latestPipelineReleases = {};

pipelines_json.remote_workflows.map(
    (pipeline) => (latestPipelineReleases[pipeline.name] = `/${pipeline.name}/${pipeline.releases[0].tag_name}/`),
);
const latestTollsURL = `/tools/docs/'+${latestToolsRelease}`;
// https://astro.build/config
export default defineConfig({
    site: 'https://nf-co.re/',
    output: 'hybrid',
    experimental: {
        assets: true,
    },
    adapter: netlify(),
    redirects: {
        [latestTollsURL]: 'https://oldsite.nf-co.re/tools/docs/latest/',
        ...latestPipelineReleases,
    },
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
        format: 'file',
    },
    vite: {
        plugins: [
            yaml(),
            FontaineTransform.vite({
                // avoid flash of unstyled text by interjecting fallback system fonts https://developer.chrome.com/blog/framework-tools-font-fallback/#using-fontaine-library
                fallbacks: ['BlinkMacSystemFont', 'Segoe UI', 'Helvetica Neue', 'Arial', 'Noto Sans'],
                resolvePath: (id) => new URL(`./public${id}`, import.meta.url),
                skipFontFaceGeneration: (fallbackName) => fallbackName === 'Font Awesome 6 Pro) fallback',
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
        shikiConfig: {
            langs: [],
            theme: githubDarkDimmed,
            wrap: false,
        },
        remarkPlugins: [emoji, remarkGfm, remarkDirective, admonitionsPlugin, mermaid, remarkMath],
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
                    const regex = /^https:\/\/(raw.)*github/;
                    if (!regex.test(url.href) && url.href?.endsWith('.md')) {
                        url.href = url.href.replace(/\.md$/, '/');
                        url.pathname = url.pathname.replace(/\.md$/, '/');
                        url.path = url.path.replace(/\.md$/, '/');
                    } else if (!regex.test(url.href) && url.href?.endsWith('.mdx')) {
                        url.href = url.href.replace(/\.mdx$/, '/');
                        url.pathname = url.pathname.replace(/\.mdx$/, '/');
                        url.path = url.path.replace(/\.mdx$/, '/');
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
                        light: 'github-light',
                    },
                },
            ],
            rehypeKatex,
        ],
    },
});
