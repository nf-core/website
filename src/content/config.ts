import { z, defineCollection } from 'astro:content';

const events = defineCollection({
    schema: z.object({
        title: z.string(),
        subtitle: z.string(),
        type: z.enum(['bytesize', 'talk', 'hackathon', 'training']),
        start_date: z.string(),
        start_time: z.string().transform((str) => str.replace(/\s+(\w+)/, ' ($1)')),
        end_date: z.string(),
        end_time: z.string().transform((str) => str.replace(/\s+(\w+)/, ' ($1)')),
        start_announcement: z.string().optional(),
        location_name: z.string().optional(),
        location_url: z.string().url().or(z.string().startsWith('#')).or(z.array(z.string().url())).optional(),
        location_latlng: z.array(z.number(), z.number()).optional(),
        address: z.string().optional(),
        start: z.date().optional(),
        end: z.date().optional(),
        duration: z.string().optional(),
        embed_at: z.string().optional(),
        import_typeform: z.boolean().optional(),
    }),
});
const docs = defineCollection({
    schema: z.object({
        title: z.string(),
        subtitle: z.string().optional(),
        weight: z.number().optional(),
        parent: z.string().optional(),
    }),
});
const about = defineCollection({
    schema: z.object({
        title: z.string(),
        description: z.string(),
        md_github_url: z.string().url().optional(),
        minHeadingDepth: z.number().optional(),
        maxHeadingDepth: z.number().optional(),
    }),
});


export const collections = {
    events: events,
    docs: docs,
    about: about,
};
