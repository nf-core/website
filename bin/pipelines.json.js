#! /usr/bin/env node
import octokit, { getDocFiles } from '../src/components/octokit.js';
import { readFileSync, writeFileSync } from 'fs';
import path from 'path';
import ProgressBar from 'progress';


// get current path
const __dirname = path.resolve();

// write the pipelines.json file
const writePipelinesJson = async () => {
  const pipelinesJson = readFileSync(path.join(__dirname, '/public/pipelines.json'));
  const namesJson = readFileSync(path.join(__dirname, '/public/pipeline_names.json'));
  const pipelines = JSON.parse(pipelinesJson);
  const names = JSON.parse(namesJson).pipeline;
  let bar = new ProgressBar('  fetching pipelines [:bar] :percent :etas', { total: names.length });
  // go through names and add or update pipelines in pipelines.json
  for (const name of names) {
    // get the details from the github repo description
    const data = await octokit.rest.repos
      .get({
        owner: 'nf-core',
        repo: name,
      })
      .then((response) => {
        // filter out entries with _url in the key name
        response.data = Object.keys(response.data)
          .filter((key) => !key.includes('_url') && !['owner', 'permissions', 'license', 'organization'].includes(key))
          .reduce((obj, key) => {
            obj[key] = response.data[key];
            return obj;
          }, {});

        return response.data;
      });

    // get number of open pull requests
    let { data: pull_requests } = await octokit.rest.pulls.list({
      owner: 'nf-core',
      repo: name,
      state: 'open',
    });
    data['open_pr_count'] = pull_requests.length;

    // get all the contributors
    let { data: contributors } = await octokit.rest.repos.listContributors({
      owner: 'nf-core',
      repo: name,
      per_page: 100,
    });

    data['contributors'] = contributors
      .filter((contrtibutor) => contrtibutor.login !== 'nf-core-bot')
      .map((contributor) => {
        return { name: contributor.login, count: contributor.contributions };
      });

    // get the releases
    let { data: releases } = await octokit.rest.repos.listReleases({
      owner: 'nf-core',
      repo: name,
    });
    // get sha of release commit
    for (const release of releases) {
      const { data: commit } = await octokit.rest.repos.getCommit({
        owner: 'nf-core',
        repo: name,
        ref: release.tag_name,
      });
      release['sha'] = commit.sha;
    }

    // get last push to dev branch
    const { data: dev_branch } = await octokit.rest.repos.listCommits({
      owner: 'nf-core',
      repo: name,
      sha: 'dev',
    });
    if (dev_branch.length > 0) {
      releases = [
        ...releases,
        { tag_name: 'dev', published_at: dev_branch[0].commit.author.date, sha: dev_branch[0].sha },
      ];
    } else {
      console.log(`No commits to dev branch found for ${name}`);
    }
    data['releases'] = releases.map(async (release) => {
      const { tag_name, published_at, sha } = release;
      const doc_files = await getDocFiles(name, release.tag_name);

      let components = await octokit
        .request('GET /repos/{owner}/{repo}/contents/{path}?ref={ref}', {
          owner: 'nf-core',
          repo: name,
          path: 'modules.json',
          ref: tag_name,
        })
        .catch((error) => {
          if (error.status === 404) {
            // console.log(`modules.json not found in ${name} ${tag_name}`);
            return;
          } else {
            console.log(error);
            return;
          }
        })
        .then((response) => {
          if (response) {
            const modules_json = JSON.parse(Buffer.from(response.data.content, 'base64').toString());
            if (modules_json.repos['nf-core/modules']) {
              if (modules_json.repos['nf-core/modules'].modules) {
                return { modules: Object.keys(modules_json.repos['nf-core/modules'].modules) };
              }
              return { modules: Object.keys(modules_json.repos['nf-core/modules']) };
            } else if (modules_json.repos['https://github.com/nf-core/modules.git']) {
              if (modules_json.repos['https://github.com/nf-core/modules.git'].subworkflows) {
                return {
                  modules: Object.keys(modules_json.repos['https://github.com/nf-core/modules.git'].modules['nf-core']),
                  subworkflows: Object.keys(
                    modules_json.repos['https://github.com/nf-core/modules.git'].subworkflows['nf-core']
                  ),
                };
              } else {
                return {
                  modules: Object.keys(modules_json.repos['https://github.com/nf-core/modules.git'].modules['nf-core']),
                };
              }
            }
          }
        });
      if (components && components.modules) {
        components.modules = components.modules.map((component) => {
          return component.replace('/', '_');
        });
      }
      return { tag_name, published_at, sha, doc_files, components };
    });

    // resolve the promises
    data['releases'] = await Promise.all(data['releases']);
    if (!pipelines.remote_workflows) {
      pipelines.remote_workflows = [];
    }
    // update in pipelines.remote_workflows if entry with name exists or add it otherwise
    const index = pipelines.remote_workflows.findIndex((workflow) => workflow.name === name);
    if (index > -1) {
      pipelines.remote_workflows[index] = data;
    } else {
      pipelines.remote_workflows.push(data);
    }
    bar.tick();
    // write the pipelines.json file
  }

  const json = JSON.stringify(pipelines, null, 4);
  await writeFileSync(path.join(__dirname, '/public/pipelines.json'), json, 'utf8');
};

writePipelinesJson();
