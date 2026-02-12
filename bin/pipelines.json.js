#! /usr/bin/env node
import { octokit, getDocFiles, getGitHubFile, githubFolderExists } from "../sites/main-site/src/components/octokit.js";
import { fetchCO2FootprintFiles } from "./s3Utils.js";

import { promises as fs, writeFileSync, existsSync } from "fs";
import yaml from "js-yaml";
import path, { join } from "path";
import ProgressBar from "progress";

// get current path
const __dirname = path.resolve();

// Helper function to extract nextflowVersion from nextflow.config content
const extractNextflowVersion = (content) => {
  if (!content) return null;
  // Look for patterns like: nextflowVersion = '!>=24.04.2'
  const pattern = /nextflowVersion\s*=\s*['"]([^'"]+)['"]/;
  const match = content.match(pattern);
  return match ? match[1] : null;
};

// Helper function to extract nf_core_version from .nf-core.yml content
const extractNfCoreVersion = (content) => {
  if (!content) return null;
  // Look for patterns like: nf_core_version: 3.3.1
  const pattern = /nf_core_version:\s*['"]?([^\s'"]+)['"]?/;
  const match = content.match(pattern);
  return match ? match[1] : null;
};

// Get pipeline name from command line argument if provided
const args = process.argv.slice(2);
const singlePipelineName = args[0] || null;

//check if pipelines.json exists
if (!existsSync(join(__dirname, "public/pipelines.json"))) {
  // create empty pipelines.json with empty remote_workflows array
  const pipelines = { remote_workflows: [] };
  const json = JSON.stringify(pipelines, null, 4);
  await writeFileSync(path.join(__dirname, "/public/pipelines.json"), json, "utf8");
}
// write the pipelines.json file
export const writePipelinesJson = async () => {
  const pipelinesJsonPromise = fs.readFile(join(__dirname, "public/pipelines.json"), "utf8");
  const ignoredTopicsPromise = fs.readFile(join(__dirname, "sites/main-site/src/config/ignored_repos.yaml"), "utf8");
  const [pipelinesJson, ignoredTopicsYaml] = await Promise.all([pipelinesJsonPromise, ignoredTopicsPromise]);

  const pipelines = JSON.parse(pipelinesJson);

  // get ignored_repos from ignored_reops.yml
  const ignored_repos = yaml.load(ignoredTopicsYaml).ignore_repos;

  // get all repos for the nf-core org or use the specified pipeline
  let names = [];
  let active_names = [];

  if (singlePipelineName) {
    names = [singlePipelineName];
    // Check if the pipeline exists
    try {
      await octokit.rest.repos.get({
        owner: "nf-core",
        repo: singlePipelineName,
      });
      active_names = [singlePipelineName];
    } catch (err) {
      console.error(`Pipeline ${singlePipelineName} not found in nf-core organization`);
      process.exit(1);
    }
    console.log(`Processing single pipeline: ${singlePipelineName}`);
  } else {
    await octokit
      .paginate(octokit.rest.repos.listForOrg, {
        org: "nf-core",
        type: "public",
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
  }

  // Only update pipeline_names.json if processing all pipelines
  if (!singlePipelineName) {
    // write pipeline_names.json
    await fs.writeFile(
      join(__dirname, "/public/pipeline_names.json"),
      JSON.stringify({ pipeline: active_names }, null),
      "utf8",
    );
  }

  // get ignored_topics from ignored_repos.yml
  const ignored_topics = yaml.load(ignoredTopicsYaml).ignore_topics;

  // Get latest tools release
  const latest_tools_release_date = (
    await octokit.rest.repos.getLatestRelease({
      owner: "nf-core",
      repo: "tools",
    })
  )?.data?.created_at;

  // Use active pipelines count for progress bar unless processing a single pipeline
  const totalPipelines = singlePipelineName ? names.length : active_names.length;
  let bar = new ProgressBar("  fetching pipelines [:bar] :percent :etas", { total: totalPipelines });

  // Process pipelines with controlled concurrency to avoid rate limits
  const CONCURRENCY_LIMIT = 5; // Process 5 pipelines concurrently
  const processPipeline = async (name) => {
    // get the details from the github repo description
    const data = await octokit.rest.repos
      .get({
        owner: "nf-core",
        repo: name,
      })
      .then((response) => {
        // filter out entries with _url in the key name
        response.data = Object.keys(response.data)
          .filter((key) => !key.includes("_url") && !["owner", "permissions", "license", "organization"].includes(key))
          .reduce((obj, key) => {
            obj[key] = response.data[key];
            return obj;
          }, {});

        return response.data;
      });

    data["allow_merge_commit"] = data.allow_merge_commit ?? -1;
    data["allow_rebase_merge"] = data.allow_rebase_merge ?? -1;
    data["allow_squash_merge"] = data.allow_squash_merge ?? -1;

    // Get team permissions
    for (const team of ["contributors", "core"]) {
      try {
        const team_permission = await octokit.request("GET /orgs/{org}/teams/{team_slug}/repos/{owner}/{repo}", {
          org: "nf-core",
          team_slug: team,
          owner: "nf-core",
          repo: name,
          headers: {
            "X-GitHub-Api-Version": "2022-11-28",
            accept: "application/vnd.github.v3.repository+json",
          },
        });
        data[`team_${team}_permission_push`] = team_permission?.data.permissions?.push ?? -1;
        data[`team_${team}_permission_admin`] = team_permission?.data.permissions?.admin ?? -1;
      } catch (err) {
        console.warn(`Failed to fetch ${team} team permission`, err);
      }
    }

    // get all rulesets for the repo with detailed info in parallel
    let ruleSetData = [];
    try {
      const ruleSets = await octokit.request("GET /repos/{owner}/{repo}/rulesets", {
        owner: "nf-core",
        repo: name,
      });

      // Fetch all ruleset details in parallel instead of sequentially
      if (ruleSets.data.length > 0) {
        ruleSetData = await Promise.all(
          ruleSets.data.map((ruleSet) =>
            octokit.request("GET /repos/{owner}/{repo}/rulesets/{ruleset_id}", {
              owner: "nf-core",
              repo: name,
              ruleset_id: ruleSet.id,
            }),
          ),
        );
      }
    } catch (err) {
      console.warn(`Failed to fetch rulesets for ${name}`, err);
    }

    // Get all branches in one API call instead of checking each individually
    let allBranches = [];
    try {
      const { data: branchesData } = await octokit.rest.repos.listBranches({
        owner: "nf-core",
        repo: name,
        per_page: 100,
      });
      allBranches = branchesData.map((branch) => branch.name);
    } catch (err) {
      console.warn(`Failed to fetch branches for ${name}`, err);
    }

    // Get branch existence & protection rules
    for (const branch of ["master", "main", "dev", "TEMPLATE"]) {
      // Check if branch exists from the list we fetched
      const branch_exists = allBranches.includes(branch);
      data[`${branch}_branch_exists`] = branch_exists;
      // No need to check rules if branch doesn't exist or is TEMPLATE (which uses an org ruleSet)
      if (!branch_exists || branch === "TEMPLATE") {
        continue;
      }

      // Initialize ALL branch protection fields with default values
      data[`${branch}_branch_protection_up_to_date`] = -1;
      data[`${branch}_branch_protection_status_checks`] = -1;
      data[`${branch}_branch_protection_required_reviews`] = -1;
      data[`${branch}_branch_protection_require_codeowner_review`] = -1;
      data[`${branch}_branch_protection_require_non_stale_review`] = -1;
      data[`${branch}_branch_protection_enforce_admins`] = -1;

      // Check if there is a ruleSet for the branch
      const ruleSet = ruleSetData.find((r) => {
        // Check if any of the include patterns match this branch
        const includePatterns = r.data.conditions?.ref_name?.include || [];
        const branchRef = "refs/heads/" + branch;
        const matches = includePatterns.includes(branchRef);
        // console.log(`Checking ruleset for ${branch}: patterns=${JSON.stringify(includePatterns)}, looking for=${branchRef}, matches=${matches}`);
        return matches;
      })?.data;
      console.log(`${name} ${branch} ruleSet:`, ruleSet ? "found" : "not found");
      if (ruleSet) {
        console.log(`Using ruleset for ${branch} branch protection`);
        data[`${branch}_uses_ruleset`] = true;

          // Check for required status checks
          const required_status_checks = ruleSet.rules.find((rule) => rule.type === "required_status_checks")
            ?.parameters?.required_status_checks;
          // console.log(`${branch} required_status_checks:`, required_status_checks);
          data[`${branch}_branch_protection_status_checks`] = required_status_checks
            ? required_status_checks.map((check) => check.context)
            : -1;

          // Check for pull request rules
          const pull_request_rule = ruleSet.rules.find((rule) => rule.type === "pull_request")?.parameters;
          // console.log(`${branch} pull_request_rule:`, pull_request_rule);
          data[`${branch}_branch_protection_required_reviews`] =
            pull_request_rule?.required_approving_review_count ?? -1;
          data[`${branch}_branch_protection_require_codeowner_review`] =
            pull_request_rule?.require_code_owner_review ?? -1;
          data[`${branch}_branch_protection_require_non_stale_review`] =
            pull_request_rule?.dismiss_stale_reviews_on_push ?? -1;

          // Check for admin enforcement
          const enforce_admins_rule = ruleSet.rules.find((rule) => rule.type === "enforce_admins");
          data[`${branch}_branch_protection_enforce_admins`] = enforce_admins_rule?.parameters?.enabled ?? -1;

          // Check that branches don't need to be up to date
          const branch_protection_up_to_date =
            ruleSet.rules.find((rule) => rule.type === "branch_protection_up_to_date")?.parameters?.enabled ?? -1;
          data[`${branch}_branch_protection_up_to_date`] = branch_protection_up_to_date ? true : false;

          console.log(`Completed ruleset processing for ${branch} branch - found ${ruleSet.rules.length} rules`);
      } else {
        // No ruleset found, check old-style branch protection rules
          // Get branch protection rules
          try {
            const branchRules = await octokit.rest.repos.getBranchProtection({
              owner: "nf-core",
              repo: name,
              branch: branch,
            });
            // console.log(`${branch} branch protection:`, branchRules?.data);
            data[`${branch}_branch_protection_status_checks`] =
              branchRules?.data?.required_status_checks?.contexts ?? -1;
            data[`${branch}_branch_protection_required_reviews`] =
              branchRules?.data?.required_pull_request_reviews?.required_approving_review_count ?? -1;
            data[`${branch}_branch_protection_require_codeowner_review`] =
              branchRules?.data?.required_pull_request_reviews?.require_code_owner_reviews ?? -1;
            data[`${branch}_branch_protection_require_non_stale_review`] =
              branchRules?.data?.required_pull_request_reviews?.dismiss_stale_reviews ?? -1;
            data[`${branch}_branch_protection_enforce_admins`] = branchRules?.data?.enforce_admins?.enabled ?? -1;
          } catch (err) {
            // Branch protection might not be configured, which is normal
            if (err.status === 404) {
              console.log(`No branch protection configured for ${branch} branch in ${name}`);
            } else {
              console.log(
                `Failed to fetch ${branch} branch protection for ${name}:`,
                err.response?.data?.message || err.message,
              );
            }
          }

      }
    }
    // remove ignored topics
    data["has_required_topics"] = ["nf-core", "nextflow", "workflow", "pipeline"].every((topic) =>
      data["topics"].includes(topic),
    );
    data["topics"] = data["topics"].filter((topic) => !ignored_topics.includes(topic));
    if (!data["has_required_topics"]) {
      data["missing_topics"] = ["nf-core", "nextflow", "workflow", "pipeline"].filter(
        (topic) => !data["topics"].includes(topic),
      );
    }

    // get number of open pull requests
    let { data: pull_requests } = await octokit.rest.pulls.list({
      owner: "nf-core",
      repo: name,
      state: "open",
    });
    data["open_pr_count"] = pull_requests.length;

    data["repository_url"] = `https://github.com/nf-core/${name}`;

    // get all the contributors
    let { data: contributors } = await octokit.rest.repos.listContributors({
      owner: "nf-core",
      repo: name,
      per_page: 100,
    });
    data["contributors"] = contributors
      .filter((contrtibutor) => contrtibutor.login !== "nf-core-bot")
      .map((contributor) => {
        return { name: contributor.login, count: contributor.contributions, avatar_url: contributor.avatar_url };
      });

    // get the releases
    let { data: releases } = await octokit.rest.repos.listReleases({
      owner: "nf-core",
      repo: name,
    });

    // remove empty values from releases (usually from draft releases)
    releases = releases.filter((release) => release.tag_name !== "" && release.draft === false);

    // remove releases that are already in the pipelines.json file
    const pipelineIndex = pipelines.remote_workflows.findIndex((workflow) => workflow.name === name);
    let new_releases = releases;

    data["released_after_tools"] =
      new Date(latest_tools_release_date.valueOf()) < new Date(releases[0]?.published_at).valueOf();

    let old_releases = [];
    if (pipelineIndex > -1) {
      old_releases = pipelines.remote_workflows[pipelineIndex].releases.filter((release) => release.tag_name !== "dev");
      const existing_releases = old_releases.map((release) => release.tag_name);
      new_releases = new_releases.filter((release) => !existing_releases.includes(release.tag_name));
    }
    // get sha and schema info for each release in parallel
    await Promise.all(
      new_releases.map(async (release) => {
        // Fetch commit info and schema existence in parallel
        const [commitData, schemaExists] = await Promise.all([
          octokit.rest.repos.getCommit({
            owner: "nf-core",
            repo: name,
            ref: release.tag_name,
          }),
          getGitHubFile(name, "nextflow_schema.json", release.tag_name).then((response) => {
            return response ? true : false;
          }),
        ]);

        release.tag_sha = commitData.data.sha;
        release.has_schema = schemaExists;
      }),
    );

    // Parallelize commit fetching and file existence checks
    const [{ data: default_branch }, lastReleaseCommit, { data: dev_branch }, isDSL2, hasNfTest, hasNfTestDev] =
      await Promise.all([
        // Get default branch commits
        octokit.rest.repos.listCommits({
          owner: "nf-core",
          repo: name,
          sha: data.default_branch,
          per_page: 1,
        }),
        // Get last release commit (only if releases exist)
        releases[0]?.tag_name
          ? octokit.rest.repos.getCommit({
              owner: "nf-core",
              repo: name,
              ref: releases[0].tag_name,
            })
          : Promise.resolve(null),
        // Get dev branch commits
        octokit.rest.repos
          .listCommits({
            owner: "nf-core",
            repo: name,
            sha: "dev",
          })
          .catch(() => ({ data: [] })), // Handle case where dev branch doesn't exist
        // Check DSL2 status
        githubFolderExists(name, "modules", releases?.[0]?.tag_name ?? data.default_branch),
        // Check nf-test usage
        getGitHubFile(name, "nf-test.config", releases?.[0]?.tag_name ?? data.default_branch).then((response) => {
          return response ? true : false;
        }),
        // Check nf-test dev usage
        getGitHubFile(name, "nf-test.config", "dev").then((response) => {
          return response ? true : false;
        }),
      ]);

    data["head_sha"] = default_branch[0]?.sha;
    data["last_release_is_head"] = lastReleaseCommit ? data["head_sha"] === lastReleaseCommit.data?.sha : false;
    data["last_release_vs_default_compare_url"] = lastReleaseCommit
      ? data["repository_url"] + "/compare/" + data["head_sha"] + "..." + lastReleaseCommit.data?.sha
      : "";

    data["commits_to_dev"] = dev_branch?.length;
    data["is_DSL2"] = isDSL2;
    data["has_nf_test"] = hasNfTest;
    data["has_nf_test_dev"] = hasNfTestDev;

    if (dev_branch.length > 0) {
      // Get dev schema in parallel with other operations
      const devSchema = await getGitHubFile(name, "nextflow_schema.json", "dev").then((response) => {
        return response ? true : false;
      });

      new_releases = [
        ...new_releases,
        {
          tag_name: "dev",
          published_at: dev_branch[0].commit.author.date,
          tag_sha: dev_branch[0].sha,
          has_schema: devSchema,
        },
      ];
    } else {
      console.log(`No commits to dev branch found for ${name}`);
    }

    // Get plugins in nextflow.config file, in the form of:
    //  plugins {
    // id 'nf-validation' // Validation of pipeline parameters and creation of an input channel from a sample sheet
    // }
    for (const branch of ["master", "main", "dev"]) {
      if (!data[`${branch}_branch_exists`]) {
        continue;
      }
      const nextflowConfig = await getGitHubFile(name, "nextflow.config", branch).then((response) => {
        if (response) {
          // Parse plugins
          const plugins = response.match(/plugins\s*{([^}]*)}/s);
          const pluginsList = plugins
            ? plugins[1]
                .split("\n")
                .filter((line) => line.includes("id"))
                .map((line) => line.match(/id\s*['"]([^'"]*)['"]/)[1])
            : [];

          // get default branch from manifest
          const manifest = {};
          const defaultBranchMatch = response.match(/manifest\s*{[^}]*defaultBranch\s*=\s*['"]([^'"]+)['"]/s);
          if (defaultBranchMatch) {
            manifest.defaultBranch = defaultBranchMatch[1];
          }

          return { plugins: pluginsList, manifest };
        }
        return { plugins: [], manifest: {} };
      });

      data[`${branch}_nextflow_config_plugins`] = nextflowConfig.plugins;
      data[`${branch}_nextflow_config_manifest`] = nextflowConfig.manifest;
    }

    new_releases = await Promise.all(
      new_releases.map(async (release) => {
        const { tag_name, published_at, tag_sha, has_schema } = release;
        const doc_files = await getDocFiles(name, release.tag_name);

        // Fetch version information
        let nextflow_version = null;
        let nf_core_version = null;

        // Fetch nextflow.config and extract nextflowVersion
        const nextflowConfig = await getGitHubFile(name, "nextflow.config", tag_name);
        if (nextflowConfig) {
          nextflow_version = extractNextflowVersion(nextflowConfig);
        }

        // Fetch .nf-core.yml and extract nf_core_version
        const nfCoreYml = await getGitHubFile(name, ".nf-core.yml", tag_name);
        if (nfCoreYml) {
          nf_core_version = extractNfCoreVersion(nfCoreYml);
        }

        // Check if this release has plugins already stored either as dev or as "main"/"master"
        let plugins = data[`main_nextflow_config_plugins`] || data[`master_nextflow_config_plugins`] || [];

        // Fetch CO2 footprint files if nf-co2footprint plugin is used and we are not in the dev tree
        console.log(plugins);
        let co2footprint_files = null;
        if (plugins.join(" ").includes("nf-co2footprint") && tag_name !== "dev") {
          console.log(`Fetching CO2 footprint files for ${name} ${tag_name}`);
          co2footprint_files = await fetchCO2FootprintFiles(name, tag_sha);
        }

        let components = await octokit
          .request("GET /repos/{owner}/{repo}/contents/{path}?ref={ref}", {
            owner: "nf-core",
            repo: name,
            path: "modules.json",
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
              const modules_json = JSON.parse(Buffer.from(response.data.content, "base64").toString());
              if (modules_json.repos["nf-core/modules"]) {
                if (modules_json.repos["nf-core/modules"].modules) {
                  return { modules: Object.keys(modules_json.repos["nf-core/modules"].modules) };
                }
                return { modules: Object.keys(modules_json.repos["nf-core/modules"]) };
              } else if (
                modules_json.repos["https://github.com/nf-core/modules.git"] &&
                modules_json.repos["https://github.com/nf-core/modules.git"].modules &&
                modules_json.repos["https://github.com/nf-core/modules.git"].modules["nf-core"]
              ) {
                if (
                  modules_json.repos["https://github.com/nf-core/modules.git"].subworkflows &&
                  modules_json.repos["https://github.com/nf-core/modules.git"].subworkflows["nf-core"]
                ) {
                  return {
                    modules: Object.keys(
                      modules_json.repos["https://github.com/nf-core/modules.git"].modules["nf-core"],
                    ),
                    subworkflows: Object.keys(
                      modules_json.repos["https://github.com/nf-core/modules.git"].subworkflows["nf-core"],
                    ),
                  };
                } else {
                  return {
                    modules: Object.keys(
                      modules_json.repos["https://github.com/nf-core/modules.git"].modules["nf-core"],
                    ),
                  };
                }
              }
            }
          });
        if (components && components.modules) {
          components.modules = components.modules.map((component) => {
            return component.replace("/", "_");
          });
        }
        return {
          tag_name,
          published_at,
          tag_sha,
          has_schema,
          doc_files,
          components,
          nextflow_version,
          nf_core_version,
          plugins,
          co2footprint_files,
        };
      }),
    );

    // Assign new_releases to data.releases
    data["releases"] = [...old_releases, ...new_releases];

    // Resolve the promises
    data["releases"] = await Promise.all(data["releases"]);

    // sort the releases by date except for dev, which should always be last
    data.releases.sort((a, b) => {
      if (a.tag_name === "dev") {
        return 1;
      } else if (b.tag_name === "dev") {
        return -1;
      } else {
        return new Date(b.published_at) - new Date(a.published_at);
      }
    });
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
    return data;
  };

  // Process pipelines with controlled concurrency
  const processInBatches = async (items, batchSize) => {
    const results = [];
    for (let i = 0; i < items.length; i += batchSize) {
      const batch = items.slice(i, i + batchSize);
      const batchResults = await Promise.all(batch.map(processPipeline));
      results.push(...batchResults);
    }
    return results;
  };

  // Process only active (non-archived) pipelines with controlled concurrency
  const pipelinesToProcess = singlePipelineName ? names.flat() : active_names.flat();
  await processInBatches(pipelinesToProcess, CONCURRENCY_LIMIT);

  // sort the pipelines by name
  pipelines.remote_workflows.sort((a, b) => {
    return a.name.localeCompare(b.name);
  });
  const json = JSON.stringify(pipelines, null, 4);

  console.log("Writing pipelines.json");
  await writeFileSync(path.join(__dirname, "/public/pipelines.json"), json, "utf8");
};

writePipelinesJson();
