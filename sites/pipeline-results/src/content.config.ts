import { defineCollection, z } from 'astro:content';
import { awsResultsLoader } from '@utils/loaders';
import pipelines_json from "@public/pipelines.json";

// Define the type for pipelines_json to match the loader requirements
type PipelineJson = {
  remote_workflows: {
    name: string;
    releases: { tag_name: string; doc_files: string[]; has_schema: boolean; tag_sha: string }[];
    archived: boolean;
    description: string;
    topics: string[];
    stargazers_count: number;
  }[];
};

// Cast the imported JSON to the expected type
const typedPipelinesJson = pipelines_json as unknown as PipelineJson;

// Define the schema for the pipeline-results collection
const pipelineResults = defineCollection({
  schema: z.object({
    pipeline: z.string(),
    version: z.string(),
    results_path: z.string(),
    content: z.array(
      z.object({
        Key: z.string(),
        LastModified: z.date(),
        ETag: z.string(),
        Size: z.number(),
        StorageClass: z.string()
      })
    ),
    commonPrefixes: z.array(
      z.object({
        Prefix: z.string()
      })
    )
  }),
  // Use the awsResultsLoader with the properly typed pipelines_json
  loader: awsResultsLoader(typedPipelinesJson)
});

// Export the collections
export const collections = {
  'pipeline-results': pipelineResults,
};

