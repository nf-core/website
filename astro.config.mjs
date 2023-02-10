// https://astro.build/config
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
import emoji from 'remark-emoji';
import remarkGfm from 'remark-gfm';
import { BUNDLED_LANGUAGES } from 'shiki';

BUNDLED_LANGUAGES = BUNDLED_LANGUAGES.map((lang) => {
    if (lang.id === 'groovy') {
        lang.aliases = ['nextflow', 'nf'];
    }
    return lang;
});

export default defineConfig({
    site: 'https://nf-co.re/',
    integrations: [svelte(), sitemap(), markdownIntegration()],
    vite: {
        plugins: [yaml()],
        ssr: {
            noExternal: ['@popperjs/core'],
        },
    },
    markdown: {
        syntaxHighlight: false,
        remarkPlugins: [emoji, remarkGfm],
        rehypePlugins: [
            rehypeSlug,
            [
                rehypeAutolinkHeadings,
                {
                    behavior: 'append',
                    content: h('i.ms-1.fas.fa-link.invisible'),
                },
            ],
            [addClasses, { table: 'table table-hover table-sm small' }],
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
        ],
    },
});