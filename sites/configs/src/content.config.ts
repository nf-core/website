import { defineCollection } from "astro:content";
import { z } from "astro/zod";
import { configLoader } from "@utils/loaders";

const configs = defineCollection({
    loader: configLoader(),
    schema: z.object({
        extension: z.literal("md").or(z.literal("config")),
        sha: z.string(),
        lastCommit: z.string().datetime(),
    }),
});

export const collections = {
    configs,
};
