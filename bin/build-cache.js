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
      const version = release.tag_name;
      for (const f of release.doc_files) {
        const cache_key = `${name}/${version}/${f}`;
        const is_cached = cache.getSync(cache_key, false) && cache.getSync(cache_key, false).length > 0;
        if (!is_cached || force) {
          await octokit
            .request('GET /repos/{owner}/{repo}/contents/{path}?ref={ref}', {
              owner: 'nf-core',
              repo: name,
              path: f,
              ref: version,
            })
            .catch((error) => {
              if (error.status === 404) {
                console.log(`File ${f} not found in ${name} ${version}`);
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
              let content = Buffer.from(response.data.content, 'base64').toString('utf-8');

              if (f.endsWith('.md')) {
                const parent_directory = f.split('/').slice(0, -1).join('/');
                // add github url to image links in markdown
                content = content.replaceAll(/!\[(.*?)\]\((.*?)\)/g, (match, p1, p2) => {
                  return `![${p1}](https://raw.githubusercontent.com/nf-core/${name}/${version}/${parent_directory}/${p2})`;
                });
                // add github url to html img tags in markdown
                content = content.replaceAll(/<img(.*?)src="(.*?)"/g, (match, p1, p2) => {
                  return `<img${p1}src="https://raw.githubusercontent.com/nf-core/${name}/${version}/${parent_directory}/${p2}"`;
                });
                // remove github warning and everything before from docs
                content = content.replace(/(.*?)(## :warning:)(.*?)(f)/s, '');
                // remove blockquote ending in "files._" from the start of the document
                content = content.replace(/(.*?)(files\._)/s, '');
                // cleanup heading
                content = content.replace('# nf-core/' + name + ': ', '# ');
                // remove everything before introduction
                content = content.replace(/.*?## Introduction/gs, '## Introduction');
              }
              cache.set(cache_key, content);
            });
        } else {
          console.log(`Already cached ${cache_key}`);
        }
      }
    }
  }
};
buildCache();
