import type { Loader } from 'astro/loaders';
import pipelines_json from '@public/pipelines.json';
import { render } from 'astro:content';

export const pipelinesLoader: Loader = {
    name: 'pipelines',
    load: async (context) => {
        pipelines_json.remote_workflows.map(async (pipeline) => {
            const readmeUrl = `https://raw.githubusercontent.com/${pipeline.full_name}/dev/README.md`;
            const response = await fetch(readmeUrl);
            const content = await response.text();

            context.store.set({
                id: pipeline.name + '-dev',
                data: {
                    name: pipeline.name,
                    content,
                },
            });
        });
    },
};
