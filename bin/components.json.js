#! /usr/bin/env node
import { octokit, getCurrentRateLimitRemaining } from '../sites/main-site/src/components/octokit.js';
import { readFileSync, writeFileSync, existsSync } from 'fs';
import path, { join } from 'path';
import ProgressBar from 'progress';
import { parse } from 'yaml';


// get current path
const __dirname = path.resolve();

// Get pipeline name from command line argument if provided
const args = process.argv.slice(2);
const singlePipelineName = args[0] || null;

console.log(await getCurrentRateLimitRemaining());
// write the components.json file
export const writeComponentsJson = async () => {
  //check if components.json exists
  if (!existsSync(join(__dirname, '/public/components.json'))) {
    // create empty components.json with empty modules and subworkflos array
    const components = { modules: [], subworkflows: [] };
    const json = JSON.stringify(components, null, 4);
    await writeFileSync(path.join(__dirname, '/public/components.json'), json, 'utf8');
  }
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
    .then((response) => response.data.tree);

  let modules = tree
    .filter((file) => file.path.includes('meta.yml') && !file.path.includes('subworkflows/'))
    .map((file) => ({
      name: file.path.replace('modules/nf-core/', '').replace('/meta.yml', '').replace('/', '_'),
      path: file.path,
      type: 'module',
    }));

  let bar = new ProgressBar('  fetching module meta.ymls [:bar] :percent :etas', { total: modules.length });

  // Fetch content for modules concurrently
  for (const module of modules) {
    const content = await octokit
      .request('GET /repos/{owner}/{repo}/contents/{path}', {
        owner: 'nf-core',
        repo: 'modules',
        path: module.path,
      })
      .then((response) => parse(Buffer.from(response.data.content, 'base64').toString()));

    module['meta'] = content;

    const index = components.modules.findIndex((m) => m.name === module.name);
    if (index > -1) {
      components.modules[index] = module;
    } else {
      components.modules.push(module);
    }
    bar.tick();
  }

  // Fetch subworkflows concurrently
  const subworkflows = tree
    .filter(
      (file) =>
        file.path.includes('meta.yml') && file.path.includes('subworkflows/') && !file.path.includes('homer/groseq'),
    )
    .map((file) => ({
      name: file.path.replace('subworkflows/nf-core/', '').replace('/meta.yml', ''),
      path: file.path,
      type: 'subworkflow',
    }));

  bar = new ProgressBar('  fetching subworkflow meta.ymls [:bar] :percent :etas', { total: subworkflows.length });

  for (const subworkflow of subworkflows) {
    const content = await octokit
      .request('GET /repos/{owner}/{repo}/contents/{path}', {
        owner: 'nf-core',
        repo: 'modules',
        path: subworkflow.path,
      })
      .then((response) => parse(Buffer.from(response.data.content, 'base64').toString()));

    subworkflow['meta'] = content;

    if (!components.subworkflows) {
      components.subworkflows = [];
    }
    const index = components.subworkflows.findIndex((m) => m.name === subworkflow.name);
    if (index > -1) {
      components.subworkflows[index] = subworkflow;
    } else {
      components.subworkflows.push(subworkflow);
    }

    if (content.modules) {
      for (const module of content.modules) {
        const index = components.modules.findIndex((m) => m.name === module);
        if (index > -1) {
          const entry = subworkflow.name;
          if (components.modules[index].subworkflows) {
            components.modules[index].subworkflows.push(entry);
          } else {
            components.modules[index].subworkflows = [entry];
          }
        }
      }
    }

    bar.tick();
  }
  // Update pipelines that use modules and subworkflows
  let pipelinesToProcess = pipelines.remote_workflows;

  if (singlePipelineName) {
    pipelinesToProcess = pipelines.remote_workflows.filter(pipeline => pipeline.name === singlePipelineName);
    if (pipelinesToProcess.length === 0) {
      console.error(`Pipeline ${singlePipelineName} not found in pipelines.json`);
      process.exit(1);
    }
    console.log(`Processing components for single pipeline: ${singlePipelineName}`);

    // Clean up existing pipeline references for this specific pipeline
    components.modules.forEach(module => {
      if (module.pipelines) {
        module.pipelines = module.pipelines.filter(p => p.name !== singlePipelineName);
        if (module.pipelines.length === 0) {
          delete module.pipelines;
        }
      }
    });

    components.subworkflows.forEach(subworkflow => {
      if (subworkflow.pipelines) {
        subworkflow.pipelines = subworkflow.pipelines.filter(p => p.name !== singlePipelineName);
        if (subworkflow.pipelines.length === 0) {
          delete subworkflow.pipelines;
        }
      }
    });
  }

  for (const pipeline of pipelinesToProcess) {
    const release = pipeline.releases[0];

    if (release.components && release.components.modules) {
      await Promise.all(
        release.components.modules.map(async (module) => {
          const index = components.modules.findIndex((m) => m.name === module);
          if (index > -1) {
            const entry = { name: pipeline.name, version: release.tag_name };
            if (components.modules[index].pipelines) {
              components.modules[index].pipelines.push(entry);
            } else {
              components.modules[index].pipelines = [entry];
            }
          }
        }),
      );
    }

    if (release.components && release.components.subworkflows) {
      await Promise.all(
        release.components.subworkflows.map(async (subworkflow) => {
          const index = components.subworkflows.findIndex((m) => m.name === subworkflow);
          if (index > -1) {
            const entry = { name: pipeline.name, version: release.tag_name };
            if (components.subworkflows[index].pipelines) {
              // check if the pipeline is already in the array
              if (!components.subworkflows[index].pipelines.some(e => e.name === entry.name && e.version === entry.version)) {
                components.subworkflows[index].pipelines.push(entry);
              }
            } else {
              components.subworkflows[index].pipelines = [entry];
            }
          }
        }),
      );
    }
  }
  // sort the modules and subworkflows by name
  components.modules.sort((a, b) => {
    if (a.name < b.name) {
      return -1;
    } else {
      return 1;
    }
  });
  components.subworkflows.sort((a, b) => {
    if (a.name < b.name) {
      return -1;
    } else {
      return 1;
    }
  });
  console.log('  writing components.json');
  // write the components.json file
  writeFileSync(path.join(__dirname, '/public/components.json'), JSON.stringify(components, null, 2));
};

writeComponentsJson();
