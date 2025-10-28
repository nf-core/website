import { z, defineCollection } from "astro:content";
import { pipelineLoader, releaseLoader } from "@utils/loaders";
import { glob } from "astro/loaders";

import pipelines_json from "@public/pipelines.json";

const pipelines = defineCollection({
    loader: pipelineLoader(pipelines_json),
});
const releases = defineCollection({
    loader: releaseLoader(pipelines_json),
});

const blog = defineCollection({
    loader: glob({ pattern: "**/[^_]*.{md,mdx}", base: "../../sites/main-site/src/content/blog" }),
});

const events = defineCollection({
    loader: glob({ pattern: "**/[^_]*.{md,mdx}", base: "../../sites/main-site/src/content/events" }),
    schema: z
        .object({
            title: z.string(),
            subtitle: z.string(),
            shortTitle: z.string().optional(),
            type: z.enum(["bytesize", "talk", "hackathon", "training"]),
            startDate: z.string().refine((s) => /^(\d{4}-\d{2}-\d{2})$/.test(s), {
                message: "startDate must be in the format YYYY-MM-DD",
            }),
            // check that it contains a time offset
            startTime: z.string().refine((s) => /^(\d{2}:\d{2})([+-]\d{2}:\d{2})$/.test(s), {
                message: "startTime must be in the format HH:MM+|-HH:MM where the +/-HH:MM is the UTC offset",
            }),
            endDate: z.string().refine((s) => /^(\d{4}-\d{2}-\d{2})$/.test(s), {
                message: "endDate must be in the format YYYY-MM-DD",
            }),
            endTime: z.string().refine((s) => /^(\d{2}:\d{2})([+-]\d{2}:\d{2})$/.test(s), {
                message: "endTime must be in the format HH:MM+|-HH:MM where the +/-HH:MM is the UTC offset",
            }),
            announcement: z
                .object({
                    text: z.string().optional(),
                    start: z.date().optional(),
                    end: z.date().optional(),
                })
                .optional(),
            locations: z
                .array(
                    z.object({
                        name: z.string().optional(),
                        links: z.string().url().or(z.string().startsWith("#")).or(z.array(z.string().url())).optional(),
                        geoCoordinates: z.array(z.number(), z.number()).optional(),
                        address: z.string().optional(),
                        country: z.string().optional(),
                        city: z.string().optional(),
                    }),
                )
                .optional(),
            links: z.array(z.string().url()).optional(),
            start: z.date().optional(),
            end: z.date().optional(),
            duration: z.string().optional(),
            embedAt: z.string().optional(),
            importTypeform: z.boolean().optional(),
            hackathonProjectListModals: z.string().optional(),
            youtubeEmbed: z.array(z.string().url()).optional().or(z.string().url()).optional(),
            hideExportButton: z.boolean().optional(),
        })
        .refine((data) => {
            // create start and end date objects
            try {
                data.start = data.start ?? new Date(data.startDate + "T" + data.startTime);
                data.end = data.end ?? new Date(data.endDate + "T" + data.endTime);
            } catch (e) {
                throw new Error("startDate and startTime must be in the format YYYY-MM-DD and HH:MM+|-HH:MM");
            }
            // check that start date is before end date
            if (data.start.getTime() > data.end.getTime()) {
                throw new Error(`startDate ${data.start} must be before endDate ${data.end}`);
            }

            // check that announcement.start is before announcement.end
            if (data.announcement?.start && data.announcement.end) {
                if (data.announcement.start.getTime() > data.announcement.end.getTime()) {
                    throw new Error("announcement.start must be before announcement.end");
                }
            }
            // check that announcement.start is set if announcement.text is
            if (data.announcement?.text && !data.announcement.start && !data.announcement.end) {
                throw new Error("announcement.start and announcement.end must be set if announcement.text is");
            }
            // check that locations country is set if locations city is set
            if (data.locations?.[0]?.city && !data.locations?.[0]?.country) {
                throw new Error("locations.country must be set if locations.city is");
            }
            // Return true if the validation should pass
            return true;
        })
        .transform((data) => {
            return {
                ...data,
                start: data.start!, // Assert dates are set after refinement
                end: data.end!,
            };
        }),
});

const advisories = defineCollection({
    loader: glob({ pattern: "**/[^_]*.{md,mdx}", base: "./src/content/advisories" }),
});

export const collections = {
    pipelines,
    releases,
    events,
    blog,
    advisories,
};
