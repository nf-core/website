#! /usr/bin/env node
import { getGitHubFile } from '../src/components/octokit.js';
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
    let releases = pipeline.releases;
    for (const release of releases) {
      release.doc_files.push('README.md'); // add the README to the cache
      release.doc_files.push('nextflow_schema.json'); // add the nextflow_schema.json to the cache
      const version = release.tag_name;
      for (let f of release.doc_files) {
        const cache_key = `${name}/${version}/${f}`;
        const is_cached = cache.getSync(cache_key, false) && cache.getSync(cache_key, false).length > 0;
        if (!is_cached || force || version === 'dev') {
          const content = await getGitHubFile(name, f, version);
          cache.set(cache_key, content);
        } else {
          console.log(`Already cached ${cache_key}`);
        }
      }
    }
  }
};
buildCache();
