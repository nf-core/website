import calloutsPlugin from './bin/remark-callouts.js';
import netlify from '@astrojs/netlify/functions';
import partytown from '@astrojs/partytown';
import prefetch from '@astrojs/prefetch';
import sitemap from '@astrojs/sitemap';
import svelte from '@astrojs/svelte';
import markdownIntegration from '@astropub/md';
import yaml from '@rollup/plugin-yaml';
import { defineConfig } from 'astro/config';
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
    site: 'https://deploy-preview-1652--nf-core.netlify.app/', // TODO: switch back to 'https://nf-co.re/'
    output: 'hybrid',
    experimental: {
        hybridOutput: true,
    },
    adapter: netlify(),
    integrations: [svelte(), sitemap(), markdownIntegration(), prefetch(), partytown()],
    vite: {
        plugins: [yaml()],
        ssr: {
            noExternal: ['@popperjs/core', 'bin/cache.js'],
        },
    },
    image: {
        service: { entrypoint: 'astro/assets/services/sharp' },
    },
    markdown: {
        syntaxHighlight: false,
        remarkPlugins: [emoji, remarkGfm, remarkDirective, calloutsPlugin],
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
                },
            ],
            // [rehypeWrap, { selector: 'pre:has(code.language-bash)', wrapper: 'div.copy-code' }],
        ],
    },
});
