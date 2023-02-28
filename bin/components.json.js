#! /usr/bin/env node
import octokit from './octokit.js';
import { readFileSync, writeFileSync } from 'fs';
import path from 'path';
import ProgressBar from 'progress';
import { parse } from 'yaml';

// get current path
const __dirname = path.resolve();

// write the pipelines.json file
const writeComponentsJson = async () => {
  const componentsJson = readFileSync(path.join(__dirname, '/public/components.json'));
  const components = JSON.parse(componentsJson);

  const pipelinesJson = readFileSync(path.join(__dirname, '/public/pipelines.json'));
  const pipelines = JSON.parse(pipelinesJson);
  // get meta.yml from nf-core/modules using octokit and git trees
  const tree = await octokit
    .request('GET /repos/{owner}/{repo}/git/trees/{tree_sha}', {
      owner: 'nf-core',
      repo: 'modules',
      tree_sha: 'master',
      recursive: 'true',
    })
    .then((response) => {
      return response.data.tree;
    });
  let modules = tree
    .filter((file) => {
      return file.path.includes('meta.yml') && !file.path.includes('subworkflows/');
    })
    .map((file) => {
      return {
        name: file.path.replace('modules/nf-core/', '').replace('/meta.yml', '').replace('/', '_'),
        path: file.path,
      };
    });
  let bar = new ProgressBar('  fetching meta.ymls [:bar] :percent :etas', { total: modules.length });

  for (const module of modules) {
    const content = await octokit
      .request('GET /repos/{owner}/{repo}/contents/{path}', {
        owner: 'nf-core',
        repo: 'modules',
        path: module.path,
      })
      .then((response) => {
        const content = parse(Buffer.from(response.data.content, 'base64').toString());
        return content;
      });
    module['meta'] = content;

    // update elements in components.modules if it exists, add it otherwise
    const index = components.modules.findIndex((m) => m.name === module.name);
    if (index > -1) {
      components.modules[index] = module;
    } else {
      components.modules.push(module);
    }
    bar.tick();
  }
  // get pipelines that use this module
  for (const pipeline of pipelines.remote_workflows) {
    const release = pipeline.releases[0];
    console.log(release.modules);
    if (release.modules) {
      for (const module of release.modules) {
        const index = components.modules.findIndex((m) => m.name === module);
        if (index > -1) {
          const entry = { name: pipeline.name, version: release.tag_name };
          if (components.modules[index].pipelines) {
            components.modules[index].pipelines.push(entry);
          } else {
            components.modules[index].pipelines = [entry];
          }
        }
      }
    }
  }
  // write the components.json file
  writeFileSync(path.join(__dirname, '/public/components.json'), JSON.stringify(components, null, 2));
};

writeComponentsJson();
