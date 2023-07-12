#! /usr/bin/env node
import { getGitHubFile, getCurrentRateLimitRemaining } from '../src/components/octokit.js';
import Cache from 'file-system-cache';
import { readFileSync } from 'fs';
import path from 'path';
import ProgressBar from 'progress';

const cache = Cache.default({
  basePath: './.cache',
  ns: 'nf-core',
});

// get current path
const __dirname = path.resolve();

// check for `--force` flag
const force = process.argv.includes('--force');

console.log(await getCurrentRateLimitRemaining());
export const buildCache = async () => {
  const pipelinesJson = readFileSync(path.join(__dirname, '/public/pipelines.json'));
  const pipelines = JSON.parse(pipelinesJson);

  let bar = new ProgressBar('  caching markdown [:bar] :percent :etas', { total: pipelines.remote_workflows.length });

  // go through the releases of each pipeline and get the files which are needed for the pipeline pages
  for (const pipeline of pipelines.remote_workflows) {
    // console.log(`Caching ${pipeline.name}`);
    const { name } = pipeline;
    let releases = pipeline.releases;
    for (const release of releases) {
      // console.log(`Caching ${name} ${release.tag_name}`);
      release.doc_files.push('README.md'); // add the README to the cache
      if (release.has_schema) {
        release.doc_files.push('nextflow_schema.json'); // add the nextflow_schema.json to the cache
      }
      const version = release.tag_name;
      await Promise.all(
        release.doc_files.map(async (f) => {
          const cache_key = `${name}/${version}/${f}`;
          // console.log(`Checking ${cache_key}`);
          const is_cached = cache.getSync(cache_key, false) && cache.getSync(cache_key, false).length > 0;
          if (!is_cached || force || version === 'dev') {
            const content = await getGitHubFile(name, f, version);
            // console.log(`Caching ${cache_key}`);
            cache.set(cache_key, content);
            // console.log(`Cached ${cache_key}`);
          } else {
            console.log(`Already cached ${cache_key}`);
          }
        })
      );
    }

  bar.tick();
  }
  return true;
};
buildCache();
