import { defineCollection } from 'astro:content';
import { pipelineLoader, releaseLoader } from '@utils/loaders'
import pipelines_json from "@public/pipelines.json";

const pipelines = defineCollection({
    loader: pipelineLoader(pipelines_json),
});
const releases = defineCollection({
    loader: releaseLoader(pipelines_json),
});
export const collections = {
    pipelines,
    releases,
};
