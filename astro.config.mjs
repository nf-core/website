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
import remarkDescription from 'astro-remark-description';
import markdownIntegration from '@astropub/md';
import icon from "astro-icon";

const latestToolsRelease = await fetch('https://api.github.com/repos/nf-core/tools/releases/latest')
    .then((res) => res.json())
    .then((json) => json.tag_name);

let latestPipelineReleases = {};
pipelines_json.remote_workflows.map(
    (pipeline) => (latestPipelineReleases[pipeline.name] = `/${pipeline.name}/${pipeline.releases[0].tag_name}/`),
);
// https://astro.build/config
export default defineConfig({
    site: 'https://nf-co.re/',
    output: 'hybrid',
    adapter: netlify(),
    redirects: {
        ...latestPipelineReleases,
    },
    integrations: [
        svelte(),
        icon({
            include: {
                // only include a subset of icons
                'file-icons': ['nextflow'],
                logos: [
                    'twitter',
                    'mastodon-icon',
                    'slack-icon',
                    'aws',
                    'microsoft-azure',
                    'github-actions',
                    'youtube-icon',
                    'linkedin',
                ],
                fa: ['github'],
                'fa-brands': ['github'],
                'line-md': ['check-list-3-twotone'],
                mdi: ['aws', 'slack', 'youtube'],
                octicon: ['link-external-16', 'table-16', 'chevron-right-16'],
                'simple-icons': ['bluesky'],
            },
        }),
        sitemap(),
        prefetch(),
        partytown({
            // Adds dataLayer.push as a forwarding-event.
            config: {
                forward: ['dataLayer.push'],
            },
        }),
        mdx(),
        markdownIntegration(),
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
                skipFontFaceGeneration: (fallbackName) => fallbackName === 'Font Awesome 6 Pro fallback',
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
        domains: ['raw.githubusercontent.com', 'unsplash.com'],
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
        remarkPlugins: [
            emoji,
            remarkGfm,
            remarkDirective,
            admonitionsPlugin,
            mermaid,
            remarkMath,
            [
                remarkDescription,
                {
                    name: 'excerpt',
                    node: (node, i, parent) => {
                        // check if parent has a child that is an html comment with the text 'end of excerpt'
                        if (
                            parent?.children?.some(
                                (child) =>
                                    (child.type === 'html' && child.value === '<!-- end of excerpt -->') ||
                                    (child.type === 'mdxFlowExpression' && child?.value === '/* end of excerpt */'),
                            )
                        ) {
                            const sibling = parent?.children[i + 1];

                            return (
                                (sibling?.type === 'html' && sibling?.value === '<!-- end of excerpt -->') ||
                                (sibling?.type === 'mdxFlowExpression' && sibling?.value === '/* end of excerpt */')
                            );
                        } else {
                            // return the first paragraph otherwise

                            // get the index of the first paragraph
                            const firstParagraphIndex = parent?.children.findIndex(
                                (child) => child.type === 'paragraph',
                            );
                            // if the node is the first paragraph, return true
                            return i === firstParagraphIndex;
                        }
                    },
                },
            ],
        ],
        // NOTE: Also update the plugins in `src/components/Markdown.svelte`!
        rehypePlugins: [
            rehypeSlug,
            [
                rehypeAutolinkHeadings,
                {
                    behavior: 'append',
                    content: h('i.ms-1.fas.fa-link.fa-xs.invisible'),
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
                    defaultLang: 'plaintext',
                    keepBackground: true,
                    theme: {
                        dark: 'github-dark',
                        light: 'github-light',
                    },
                },
            ],
            rehypeKatex,
        ],
    },
});
