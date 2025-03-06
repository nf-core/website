import { defineCollection } from 'astro:content';
import { pipelineLoader, releaseLoader } from '@utils/loaders'
import pipelines_json from "@public/pipelines.json";

// Define the expected type for pipelines_json to match the loader requirements
type PipelinesJson = {
    remote_workflows: {
        name: string;
        releases: { tag_name: string; doc_files: string[]; has_schema: boolean }[];
        archived: boolean;
        description: string;
        topics: string[];
        stargazers_count: number;
    }[];
};

// Cast the imported JSON to the expected type
const typedPipelinesJson = pipelines_json as unknown as PipelinesJson;

const pipelines = defineCollection({
    loader: pipelineLoader(typedPipelinesJson),
});
const releases = defineCollection({
    loader: releaseLoader(typedPipelinesJson),
});
export const collections = {
    pipelines,
    releases,
};
