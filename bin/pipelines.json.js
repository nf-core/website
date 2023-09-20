#! /usr/bin/env node
import octokit, { getDocFiles, getCurrentRateLimitRemaining, getGitHubFile } from '../src/components/octokit.js';
import { promises as fs, writeFileSync, existsSync } from 'fs';
import yaml from 'js-yaml';
import path, { join } from 'path';
import ProgressBar from 'progress';


// get current path
const __dirname = path.resolve();

console.log(await getCurrentRateLimitRemaining());
//check if pipelines.json exists
if (!existsSync(join(__dirname, 'public/pipelines.json'))) {
  // create empty pipelines.json with empty remote_workflows array
  const pipelines = { remote_workflows: [] };
  const json = JSON.stringify(pipelines, null, 4);
  await writeFileSync(path.join(__dirname, '/public/pipelines.json'), json, 'utf8');
}
// write the pipelines.json file
export const writePipelinesJson = async () => {
  const pipelinesJsonPromise = fs.readFile(join(__dirname, 'public/pipelines.json'), 'utf8');
  const ignoredTopicsPromise = fs.readFile(join(__dirname, 'src/config/ignored_repos.yaml'), 'utf8');
  const [pipelinesJson, ignoredTopicsYaml] = await Promise.all([pipelinesJsonPromise, ignoredTopicsPromise]);

  const pipelines = JSON.parse(pipelinesJson);

  // get ignored_repos from ignored_reops.yml
  const ignored_repos = yaml.load(ignoredTopicsYaml).ignore_repos;

  // get all repos for the nf-core org
  let names = await octokit
    .paginate(octokit.rest.repos.listForOrg, {
      org: 'nf-core',
      type: 'public',
      per_page: 100,
    })
    .then((response) => {
      // filter out repos that are in the ignored_repos list
      return response
        .filter((repo) => !ignored_repos.includes(repo.name))
        .map((repo) => repo.name)
        .sort();
    });
  // write pipeline_names.json
  await fs.writeFile(
    join(__dirname, '/public/pipeline_names.json'),
    JSON.stringify({ pipeline: names }, null, 4),
    'utf8',
  );

  // get ignored_topics from ignored_reops.yml
  const ignored_topics = yaml.load(ignoredTopicsYaml).ignore_topics;

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
    // remove ignored topics
    data['topics'] = data['topics'].filter((topic) => !ignored_topics.includes(topic));
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
        return { name: contributor.login, count: contributor.contributions, avatar_url: contributor.avatar_url };
      });

    // get the releases
    let { data: releases } = await octokit.rest.repos.listReleases({
      owner: 'nf-core',
      repo: name,
    });

    // remove releases that are already in the pipelines.json file
    const index = pipelines.remote_workflows.findIndex((workflow) => workflow.name === name);
    let new_releases = releases;

    let old_releases = [];
    if (index > -1) {
      old_releases = pipelines.remote_workflows[index].releases.filter((release) => release.tag_name !== 'dev');
      const existing_releases = old_releases.map((release) => release.tag_name);
      new_releases = new_releases.filter((release) => !existing_releases.includes(release.tag_name));
    }
    // get sha for each release (needed for aws viewer)
    for (const release of new_releases) {
      const { data: commit } = await octokit.rest.repos.getCommit({
        owner: 'nf-core',
        repo: name,
        ref: release.tag_name,
      });
      release.tag_sha = commit.sha;
      // check if schema file exists for release
      release.has_schema = await getGitHubFile(name, 'nextflow_schema.json', release.tag_name).then((response) => {
        return response ? true : false;
      });
    }

    // get last push to dev branch
    const { data: dev_branch } = await octokit.rest.repos.listCommits({
      owner: 'nf-core',
      repo: name,
      sha: 'dev',
    });
    if (dev_branch.length > 0) {
      new_releases = [
        ...new_releases,
        {
          tag_name: 'dev',
          published_at: dev_branch[0].commit.author.date,
          tag_sha: dev_branch[0].sha,
          has_schema: await getGitHubFile(name, 'nextflow_schema.json', 'dev').then((response) => {
            return response ? true : false;
          }),
        },
      ];
    } else {
      console.log(`No commits to dev branch found for ${name}`);
    }
    new_releases = await Promise.all(
      new_releases.map(async (release) => {
        const { tag_name, published_at, tag_sha, has_schema } = release;
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
                if (
                  modules_json.repos['https://github.com/nf-core/modules.git'].subworkflows &&
                  modules_json.repos['https://github.com/nf-core/modules.git'].subworkflows['nf-core']
                ) {
                  return {
                    modules: Object.keys(
                      modules_json.repos['https://github.com/nf-core/modules.git'].modules['nf-core'],
                    ),
                    subworkflows: Object.keys(
                      modules_json.repos['https://github.com/nf-core/modules.git'].subworkflows['nf-core'],
                    ),
                  };
                } else {
                  return {
                    modules: Object.keys(
                      modules_json.repos['https://github.com/nf-core/modules.git'].modules['nf-core'],
                    ),
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
        return { tag_name, published_at, tag_sha, has_schema, doc_files, components };
      }),
    );

    // Assign new_releases to data.releases
    data['releases'] = [...old_releases, ...new_releases];

    // Resolve the promises
    data['releases'] = await Promise.all(data['releases']);

    // sort the releases by date except for dev, which should always be last
    data.releases.sort((a, b) => {
      if (a.tag_name === 'dev') {
        return 1;
      } else if (b.tag_name === 'dev') {
        return -1;
      } else {
        return new Date(b.published_at) - new Date(a.published_at);
      }
    });
    if (!pipelines.remote_workflows) {
      pipelines.remote_workflows = [];
    }
    // update in pipelines.remote_workflows if entry with name exists or add it otherwise
    if (index > -1) {
      pipelines.remote_workflows[index] = data;
    } else {
      pipelines.remote_workflows.push(data);
    }
    bar.tick();
  }

  // sort the pipelines by name
  pipelines.remote_workflows.sort((a, b) => {
    return a.name.localeCompare(b.name);
  });
  const json = JSON.stringify(pipelines, null, 4);

  console.log('Writing pipelines.json');
  await writeFileSync(path.join(__dirname, '/public/pipelines.json'), json, 'utf8');
};

writePipelinesJson();
