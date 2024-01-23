import { z, defineCollection } from 'astro:content';

const events = defineCollection({
    schema: z
        .object({
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
            announcement: z
                .object({
                    text: z.string().optional(),
                    start: z.date().optional(),
                    end: z.date().optional(),
                })
                .optional(),
            location_name: z.string().optional(),
            location_url: z.string().url().or(z.string().startsWith('#')).or(z.array(z.string().url())).optional(),
            location_latlng: z.array(z.number(), z.number()).optional(),
            address: z.string().optional(),
            start: z.date().optional(),
            end: z.date().optional(),
            duration: z.string().optional(),
            embed_at: z.string().optional(),
            import_typeform: z.boolean().optional(),
        })
        .refine((data) => {
            // create start and end date objects
            try {
                data.start = new Date(data.start_date + 'T' + data.start_time);
                data.end = new Date(data.end_date + 'T' + data.end_time);
            } catch (e) {
                throw new Error('start_date and start_time must be in the format YYYY-MM-DD and HH:MM+|-HH:MM');
            }
            // check that start date is before end date
            if (data.start.getTime() > data.end.getTime()) {
                throw new Error(`start_date ${data.start} must be before end_date ${data.end}`);
            }

            // check that announcement.start is before announcement.end
            if (data.announcement?.start && data.announcement.end) {
                if (data.announcement.start.getTime() > data.announcement.end.getTime()) {
                    throw new Error('announcement.start must be before announcement.end');
                }
            }
            // check that announcement.start is set if announcement.text is
            if (data.announcement?.text && !data.announcement.start && !data.announcement.end) {
                throw new Error('announcement.start and announcement.end must be set if announcement.text is');
            }
            // Return true if the validation should pass
            return true;
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
    schema: z
        .object({
            title: z.string(),
            subtitle: z.string(),
            headerImage: z.string().url().optional(),
            headerImageAlt: z.string().optional(),
            label: z.array(z.string()),
            pubDate: z.date(),
            authors: z.array(z.string()),
            announcement: z
                .object({
                    text: z.string().optional(),
                    start: z.date().optional(),
                    end: z.date().optional(),
                })
                .optional(),
        })
        .refine((data) => {
            // Check if headerImage is present but headerImageAlt is not
            if (data.headerImage && !data.headerImageAlt) {
                throw new Error('Please provide alt text for your `headerImage` in `headerImageAlt`.');
            }
            // check that announcement.start is before announcement.end
            if (data.announcement?.start && data.announcement.end) {
                if (data.announcement.start.getTime() > data.announcement.end.getTime()) {
                    throw new Error('`announcement.start` must be before `announcement.end`');
                }
            }
            // check that announcement.start is set if announcement.text is
            if (data.announcement?.text && !data.announcement.start && !data.announcement.end) {
                throw new Error('`announcement.start` and `announcement.end` must be set if `announcement.text` is');
            }
            // Return true if the validation should pass
            return true;
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
