#! /usr/bin/env node
import octokit, { getDocFiles, getCurrentRateLimitRemaining, getGitHubFile, githubFolderExists } from '../src/components/octokit.js';
import { promises as fs, writeFileSync, existsSync } from 'fs';
import yaml from 'js-yaml';
import path, { join } from 'path';
import ProgressBar from 'progress';
import cache from './cache.js';

// get current path
const __dirname = path.resolve();


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
  let names = [];
  let active_names = [];

  await octokit
    .paginate(octokit.rest.repos.listForOrg, {
      org: 'nf-core',
      type: 'public',
      per_page: 100,
    })
    .then((response) => {
      names.push(
        response
          .filter((repo) => !ignored_repos.includes(repo.name))
          .map((repo) => repo.name)
          .sort(),
      );
      active_names.push(
        response
          .filter((repo) => !ignored_repos.includes(repo.name) && !repo.archived)
          .map((repo) => repo.name)
          .sort(),
      );
    });
  names = names.flat();
  active_names = active_names.flat();
  // write pipeline_names.json
  await fs.writeFile(
    join(__dirname, '/public/pipeline_names.json'),
    JSON.stringify({ pipeline: active_names }, null),
    'utf8',
  );

  // get ignored_topics from ignored_reops.yml
  const ignored_topics = yaml.load(ignoredTopicsYaml).ignore_topics;

  // Get latest tools release
  const latest_tools_release_date = (await octokit.rest.repos.getLatestRelease({
    owner: 'nf-core',
    repo: 'tools',
  }))?.data?.created_at;

  let bar = new ProgressBar('  fetching pipelines [:bar] :percent :etas', { total: names.length });

  // go through names and add or update pipelines in pipelines.json
  for (const name of names.flat()?.splice(0,1)) {
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

    const repoInfo = await octokit.rest.repos.get({
      owner: 'nf-core',
      repo: name,
    });


    data['allow_merge_commit'] = repoInfo.data.allow_merge_commit ?? -1;
    data['allow_rebase_commit'] = repoInfo.data.allow_rebase_merge ?? -1;
    data['allow_squash_commit'] = repoInfo.data.allow_squash_merge ?? -1;

    // Get branch existence & protection rules
    for(const branch of [
      'main',
      'dev',
      'TEMPLATE'
    ]) {
      // Initialize to -1 (unknown)
      data[`${branch}_branch_exists`] = -1;
      data[`${branch}_branch_protection_up_to_date`] = -1;
      data[`${branch}_branch_protection_status_checks`] = -1;
      data[`${branch}_branch_protection_required_reviews`] = -1;
      data[`${branch}_branch_protection_require_codeowner_review`] = -1;
      data[`${branch}_branch_protection_require_non_stale_review`] = -1;
      data[`${branch}_branch_protection_enforce_admins`] = -1;
      data[`${branch}_restrict_push`] = -1;

      // Check if branch exists
      try {
        const branch_exists = await octokit.rest.repos.getBranch({
          owner: 'nf-core',
          repo: name,
          branch: 'TEMPLATE',
        }).then(() => true).catch((err) => {
          if (err.status === 404) {
            return false;
          }
          throw err;
        });
        data[`${branch}_branch_exists`] = branch_exists;
      } catch(err) {
        console.warn(`Failed to fetch ${branch} branch`, err);
      }

      if(branch !== 'TEMPLATE') {
        // Get branch protection rules
        try {
          const rules = await octokit.rest.repos.getBranchProtection({
            owner: 'nf-core',
            repo: name,
            branch: branch === 'dev' ? 'dev' : data.default_branch,
          });
          data[`${branch}_branch_protection_up_to_date`] = rules?.data?.res
          data[`${branch}_branch_protection_status_checks`] = rules?.data?.required_status_checks ?? -1;
          data[`${branch}_branch_protection_required_reviews`] = rules?.data?.required_pull_request_reviews?.required_approving_review_count ?? -1;
          data[`${branch}_branch_protection_require_codeowner_review`] = rules?.data?.required_pull_request_reviews?.require_code_owner_reviews ?? -1;
          data[`${branch}_branch_protection_require_non_stale_review`] = rules?.data?.required_pull_request_reviews?.dismiss_stale_reviews ?? -1;
          data[`${branch}_branch_protection_enforce_admins`] = rules?.data?.enforce_admins?.enabled ?? -1;
        } catch(err)  {
          console.log(`Failed to fetch ${branch} branch protection`, err);
        }
      } else {
        // Template branch protection rules
        try {
          restrictions = await octokit.rest.repos.getBranchProtection({
            owner: 'nf-core',
            repo: name,
            branch: 'TEMPLATE',
          });
          data[`${branch}_restrict_push`] = restrictions?.users?.length === 1 && restrictions?.users?.[0]?.login === 'nf-core-bot' ? true : false;
        } catch(err)  {
          console.log(`Failed to fetch ${branch} branch push restrictions`, err);
        }
      }
    }
    // remove ignored topics
    data['topics'] = data['topics'].filter((topic) => !ignored_topics.includes(topic));
    data['has_required_topics'] = ['nf-core', 'nextflow', 'workflow', 'pipeline'].every((topic) => data['topics'].includes(topic));
    // get number of open pull requests
    let { data: pull_requests } = await octokit.rest.pulls.list({
      owner: 'nf-core',
      repo: name,
      state: 'open',
    });
    data['open_pr_count'] = pull_requests.length;

    data['repository_url'] = `https://github.com/nf-core/${name}`;

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

    data['released_after_tools'] = new Date(latest_tools_release_date.valueOf()) < new Date(releases[0]?.published_at).valueOf();

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

    // get last push to main branch and compare to that of the newest release
    const { data: default_branch } = await octokit.rest.repos.listCommits({
      owner: 'nf-core',
      repo: name,
      sha: data.default_branch,
      per_page: 1,
    });
    data['head_sha'] = default_branch[0]?.sha;

    const lastReleaseCommit = await octokit.rest.repos.getCommit({
      owner: 'nf-core',
      repo: name,
      ref: releases[0]?.tag_name,
    });
    data['last_release_is_head'] = data['head_sha'] === lastReleaseCommit.data?.sha;
    data['last_release_vs_default_compare_url'] = data['repository_url'] + '/compare/' + data['head_sha'] + '...' + lastReleaseCommit.data?.sha;

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

    // Get DSL2 status
    data['is_DSL2'] = await githubFolderExists(name, 'modules', releases?.[0]?.tag_name ?? data.default_branch);

    // Get nf-test uage
    data['has_nf_test'] = await getGitHubFile(name, 'nf-test.config', releases?.[0]?.tag_name ?? data.default_branch).then((response) => {
      return response ? true : false;
    });
    data['has_nf_test_dev'] = await getGitHubFile(name, 'nf-test.config', 'dev').then((response) => {
      return response ? true : false;
    });

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

        // cache release body except for dev
        if (release.tag_name !== 'dev') {
          const cache_key = `${name}/${release.tag_name}/body`;
          // const is_cached = cache.getSync(cache_key, false) && cache.getSync(cache_key, false).length > 0;
          const is_cached = false;
          if (!is_cached) {
            // wrap github urls in markdown links if they are to the same repo and not already inside a link
            release.body = release.body.replaceAll(
              /(?<!\]\()https:\/\/github\.com\/nf-core\/([^\/]+)\/([^\/]+)\/([^\/\n]*)(?![\)\]])/g,
              (match, p1, p2, p3) => {
                if (p1 === name && ['pull', 'issues', 'compare'].includes(p2)) {
                  const prefix = p2 !== 'compare' ? '#' : '';
                  return `[${prefix}${p3}](${match})`;
                }
                return match;
              },
            );
            // replace usernames with links to github profiles
            release.body = release.body.replaceAll(/@(\w+([-]\w+)*)/g, (match, p1) => {
              return `[${match}](https://github.com/${p1})`;
            });
            cache.set(cache_key, release.body);
          }
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
