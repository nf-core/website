// https://astro.build/config
import svelte from '@astrojs/svelte';
import yaml from '@rollup/plugin-yaml';
import { defineConfig } from 'astro/config';
import emoji from 'remark-emoji';
import remarkGfm from 'remark-gfm';
import remarkToc from 'remark-toc';

// https://astro.build/config
export default defineConfig({
    experimental: {
        contentCollections: true,
        errorOverlay: true,
    },
    site: 'https://nf-co.re/',
    integrations: [svelte()],
    vite: {
        plugins: [yaml()],
        ssr: {
            noExternal: ['@popperjs/core', 'svelte-toc'],
        },
    },
    markdown: {
        remarkPlugins: [remarkToc, emoji, remarkGfm],
    },
});
