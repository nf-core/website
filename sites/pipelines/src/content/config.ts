import { pipelinesLoader } from '@main-site/loaders/pipelines.ts';
import { defineCollection } from 'astro:content';
import pipelines_json from '@public/pipelines.json';

const pipeline = defineCollection({
    loader: pipelinesLoader,
});

export const collections = {
    pipeline: pipeline,
};
