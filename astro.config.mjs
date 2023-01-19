// https://astro.build/config
import sitemap from '@astrojs/sitemap';
import svelte from '@astrojs/svelte';
import yaml from '@rollup/plugin-yaml';
import { defineConfig } from 'astro/config';
import { h } from 'hastscript';
import groovy from 'highlight.js/lib/languages/groovy';
import shell from 'highlight.js/lib/languages/shell';
import addClasses from 'rehype-add-classes';
import rehypeAutolinkHeadings from 'rehype-autolink-headings';
import rehypeHighlight from 'rehype-highlight';
import rehypeSlug from 'rehype-slug';
import emoji from 'remark-emoji';
import remarkGfm from 'remark-gfm';

export default defineConfig({
    experimental: {
        contentCollections: true,
        errorOverlay: true,
    },
    site: 'https://nf-co.re/',
    integrations: [svelte(), sitemap()],
    vite: {
        plugins: [yaml()],
        ssr: {
            noExternal: ['@popperjs/core', 'svelte-toc'],
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
                rehypeHighlight,
                { languages: { groovy, shell }, aliases: { groovy: 'nextflow', shell: 'console', shell: 'git' } },
            ],
        ],
    },
});
