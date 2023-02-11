#! /usr/bin/env node
import octokit from './octokit.js';
import Cache from 'file-system-cache';
import { readFileSync } from 'fs';
import path from 'path';

const cache = Cache.default({
  basePath: './.cache',
  ns: 'nf-core',
});

// get current path
const __dirname = path.resolve();

// check for `--force` flag
const force = process.argv.includes('--force');

const buildCache = async () => {
  const pipelinesJson = readFileSync(path.join(__dirname, '/public/pipelines.json'));
  const pipelines = JSON.parse(pipelinesJson);
  // go through the releases of each pipeline and get the files which are needed for the pipeline pages
  for (const pipeline of pipelines.remote_workflows) {
    const { name } = pipeline;
    const releases = pipeline.releases;
    for (const release of releases) {
      release.doc_files.push('README.md'); // add the README to the cache
      release.doc_files.push('nextflow_schema.json'); // add the nextflow_schema.json to the cache
      for (const f of release.doc_files) {
        const cache_key = `${name}/${release.tag_name}/${f}`;
        const is_cached = cache.getSync(cache_key, false) && cache.getSync(cache_key, false).length > 0;
        if (!is_cached || force) {
          await octokit
            .request('GET /repos/{owner}/{repo}/contents/{path}?ref={ref}', {
              owner: 'nf-core',
              repo: name,
              path: f,
              ref: release.tag_name,
            })
            .catch((error) => {
              if (error.status === 404) {
                console.log(`File ${f} not found in ${name} ${release.tag_name}`);
                console.log(error.request.url);
                return;
              } else {
                console.log(error);
                return;
              }
            })
            .then((response) => {
              if (response == null) {
                return;
              }
              console.log('Caching ', cache_key);
              cache.set(cache_key, response.data.content);
            });
        } else {
          console.log(`Already cached ${cache_key}`);
        }
      }
    }
  }
};
buildCache();
