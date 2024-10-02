import type { Loader } from 'astro/loaders';
import pipelines_json from '@public/pipelines.json';
import { render } from 'astro:content';

export const pipelinesLoader: Loader = {
    name: 'pipelines',
    load: async ({ store, logger, parseData, meta, generateDigest }) => {
        pipelines_json.remote_workflows.map(async (pipeline) => {
            const readmeUrl = `https://raw.githubusercontent.com/${pipeline.full_name}/dev/README.md`;
            const response = await fetch(readmeUrl);
            const content = await response.text();
            logger.info('===================================');
            store.clear();
            store.set({
                id: pipeline.name + '-dev2',
                digest: generateDigest(content),
                body: content,
                data: {
                    name: pipeline.name,
                },
            });
        });
    },
};
