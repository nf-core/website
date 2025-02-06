import { z, defineCollection } from 'astro:content';
import { glob } from 'astro/loaders';

const pipelines = defineCollection({
    loader: glob({ pattern: '**/[^_]*.{md,mdx}', base: './src/content/pipelines' }),
});

export const collections = {
    pipelines,
};
