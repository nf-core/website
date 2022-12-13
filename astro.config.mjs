// https://astro.build/config
import svelte from "@astrojs/svelte";
import { defineConfig } from 'astro/config';

// https://astro.build/config
export default defineConfig({
    site: 'https://nf-co.re/',
    integrations: [svelte()],
});
