import { defineCollection } from 'astro:content';
import { pipelineLoader} from '@utils/loaders'
import pipelines_json from "@public/pipelines.json";

const pipelines = defineCollection({
    loader: pipelineLoader(pipelines_json),
});

export const collections = {
    pipelines,
};
