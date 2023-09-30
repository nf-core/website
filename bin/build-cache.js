#! /usr/bin/env node
import { getGitHubFile, getCurrentRateLimitRemaining } from '../src/components/octokit.js';
import { readFileSync, writeFileSync, existsSync, mkdirSync } from 'fs';
import path from 'path';
import ProgressBar from 'progress';

// get current path
const __dirname = path.resolve();

// check for `--force` flag
const force = process.argv.includes('--force');
(async () => {
  console.log(await getCurrentRateLimitRemaining());
  const buildCache = async () => {
    // build the pipeline cache
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
        for (const f of release.doc_files) {
          const cache_key = `${name}/${version}/${f}`;
          // console.log(`Checking ${cache_key}`);
          const is_cached = existsSync(path.join(__dirname, '.cache', cache_key));
          if (!is_cached || force || version === 'dev') {
            const content = await getGitHubFile(name, f, version);
            // console.log(`Caching ${cache_key}`);
            // cache.set(cache_key, content);
            // console.log(`Cached ${cache_key}`);
            //generate folder structure
            const parent = cache_key.split('/').slice(0, -1).join('/');
            mkdirSync(path.join(__dirname, '.cache', parent), { recursive: true });
            writeFileSync(path.join(__dirname, '.cache', cache_key), content);
          } else {
            // console.log(`Already cached ${cache_key}`);
          }
        }
      }

      bar.tick();
    }

    console.log('Done');
    return true;
  };
  await buildCache();
})();
