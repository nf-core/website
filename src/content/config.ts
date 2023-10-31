import { z, defineCollection } from 'astro:content';

const events = defineCollection({
    schema: z.object({
        title: z.string(),
        subtitle: z.string(),
        type: z.enum(['bytesize', 'talk', 'hackathon', 'training']),
        start_date: z.string().refine((s) => /^(\d{4}-\d{2}-\d{2})$/.test(s), {
            message: 'start_date must be in the format YYYY-MM-DD',
        }),
        // check that it contains a time offset
        start_time: z.string().refine((s) => /^(\d{2}:\d{2})([+-]\d{2}:\d{2})$/.test(s), {
            message: 'start_time must be in the format HH:MM+|-HH:MM where the +/-HH:MM is the UTC offset',
        }),
        end_date: z.string().refine((s) => /^(\d{4}-\d{2}-\d{2})$/.test(s), {
            message: 'end_date must be in the format YYYY-MM-DD',
        }),
        end_time: z.string().refine((s) => /^(\d{2}:\d{2})([+-]\d{2}:\d{2})$/.test(s), {
            message: 'end_time must be in the format HH:MM+|-HH:MM where the +/-HH:MM is the UTC offset',
        }),
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

const blog = defineCollection({
    schema: z.object({
        title: z.string(),
        subtitle: z.string(),
        headerImage: z.string().url().optional(),
        label: z.array(z.string()),
        pubDate: z.date(),
    }),
});

const pipelines = defineCollection({});

export const collections = {
    events: events,
    docs: docs,
    about: about,
    pipelines: pipelines,
    blog: blog,
};
