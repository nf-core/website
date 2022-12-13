// https://astro.build/config

import { defineConfig } from 'astro/config';
import svelte from "@astrojs/svelte";
import yaml from '@rollup/plugin-yaml';

// https://astro.build/config
export default defineConfig({
    site: 'https://nf-co.re/',
    integrations: [svelte()],
    vite: {
        plugins: [yaml()],
    },
});
