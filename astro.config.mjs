// https://astro.build/config
import sitemap from '@astrojs/sitemap';
import svelte from '@astrojs/svelte';
import markdownIntegration from '@astropub/md';
import yaml from '@rollup/plugin-yaml';
import { defineConfig } from 'astro/config';
import { h } from 'hastscript';
import groovy from 'highlight.js/lib/languages/groovy';
import shell from 'highlight.js/lib/languages/shell';
import addClasses from 'rehype-add-classes';
import rehypeAutolinkHeadings from 'rehype-autolink-headings';
import rehypeHighlight from 'rehype-highlight';
import rehypeSlug from 'rehype-slug';
import urls from 'rehype-urls';
import emoji from 'remark-emoji';
import remarkGfm from 'remark-gfm';

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
            // [
            //     rehypeHighlight,
            //     { languages: { groovy, shell,bibtex }, aliases: { groovy: 'nextflow', shell: 'console', shell: 'git' } },
            // ],
        ],
    },
});
