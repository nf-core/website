import { z, defineCollection } from 'astro:content';
import { glob } from 'astro/loaders';

const docs = defineCollection({
    loader: glob({ pattern: '**/[^_]*.{md,mdx}', base: './src/content/docs' }),
    schema: z.object({
        title: z.string(),
        subtitle: z.string().optional(),
        shortTitle: z.string().optional(),
        weight: z.number().optional(),
        parent: z.string().optional(),
        parentWeight: z.number().optional(),
        type: z.enum(['tutorial']).optional(),
        markdownPlugin: z.enum(['checklist', 'addNumbersToHeadings']).optional(),
    }),
});

const api_reference = defineCollection({
    loader: glob({ pattern: '**/[^_]*.{md,mdx}', base: './src/content/api_reference' }),
});

export const collections = {
    docs,
    api_reference,
};
