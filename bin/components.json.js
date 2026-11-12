#! /usr/bin/env node
import { octokit, getCurrentRateLimitRemaining, getGitHubFile } from "../sites/main-site/src/components/octokit.js";
import { readFileSync, writeFileSync, existsSync } from "fs";
import path, { join } from "path";
import ProgressBar from "progress";
import { parse } from "yaml";

// get current path
const __dirname = path.resolve();

// Get pipeline name from command line argument if provided
const args = process.argv.slice(2);
const singlePipelineName = args[0] || null;

// cap concurrent GitHub requests; a cold cache would otherwise fire thousands at once
// and trip GitHub's secondary rate limits
const createLimiter = (max) => {
  let active = 0;
  const queue = [];
  const next = () => {
    active--;
    if (queue.length > 0) {
      queue.shift()();
    }
  };
  return (task) =>
    new Promise((resolve, reject) => {
      const run = () => {
        active++;
        task().then(resolve, reject).finally(next);
      };
      if (active < max) {
        run();
      } else {
        queue.push(run);
      }
    });
};
const limit = createLimiter(15);

console.log(await getCurrentRateLimitRemaining());
// write the components.json file
export const writeComponentsJson = async () => {
  //check if components.json exists
  if (!existsSync(join(__dirname, "/public/components.json"))) {
    // create empty components.json with empty modules and subworkflos array
    const components = { modules: [], subworkflows: [] };
    const json = JSON.stringify(components, null, 4);
    await writeFileSync(path.join(__dirname, "/public/components.json"), json, "utf8");
  }
  const componentsJson = readFileSync(path.join(__dirname, "/public/components.json"));
  const components = JSON.parse(componentsJson);

  const pipelinesJson = readFileSync(path.join(__dirname, "/public/pipelines.json"));
  const pipelines = JSON.parse(pipelinesJson);

  // get the commit SHA of master HEAD — pins the tree listing, raw-content fetches,
  // and the emitted commit_sha to one consistent snapshot
  const commitSha = await octokit.rest.repos
    .listCommits({
      owner: "nf-core",
      repo: "modules",
      sha: "master",
      per_page: 1,
    })
    .then((response) => response.data[0].sha);

  // get meta.yml from nf-core/modules using octokit and git trees
  const tree = await octokit
    .request("GET /repos/{owner}/{repo}/git/trees/{tree_sha}", {
      owner: "nf-core",
      repo: "modules",
      tree_sha: commitSha,
      recursive: "true",
    })
    .then((response) => response.data.tree);

  const metaYmlFiles = tree.filter((f) => /^(modules|subworkflows)\/nf-core\/.*\/meta\.yml$/.test(f.path));
  const componentRoots = new Set(metaYmlFiles.map((f) => f.path.replace("/meta.yml", "")));
  // tree SHA of each component directory — changes when *any* file in the component changes,
  // so it doubles as the cache key for meta/git_sha
  const treeShaByRoot = new Map(
    tree.filter((e) => e.type === "tree" && componentRoots.has(e.path)).map((e) => [e.path, e.sha]),
  );

  // Pre-build a map from component root directory → all its file paths (O(n) single pass)
  const filesByRoot = new Map();
  for (const entry of tree) {
    if (entry.type !== "blob") continue;
    const parts = entry.path.split("/");
    // walk up from the deepest ancestor directory until we hit a component root
    // (roots are never shallower than modules/nf-core/<tool>)
    for (let depth = parts.length - 1; depth >= 3; depth--) {
      const candidate = parts.slice(0, depth).join("/");
      if (componentRoots.has(candidate)) {
        if (!filesByRoot.has(candidate)) filesByRoot.set(candidate, []);
        filesByRoot.get(candidate).push(entry.path);
        break;
      }
    }
  }

  const buildComponent = (file, type) => {
    const root = file.path.replace("/meta.yml", "");
    const name =
      type === "module"
        ? root.replace("modules/nf-core/", "").replace("/", "_")
        : root.replace("subworkflows/nf-core/", "");
    return {
      name,
      path: file.path,
      tree_sha: treeShaByRoot.get(root),
      files: filesByRoot.get(root) ?? [],
      type,
    };
  };

  // fill in meta (parsed meta.yml) and git_sha (last commit touching the component),
  // reusing the existing entry when the component directory is unchanged
  const enrichComponent = (component, existing, bar) =>
    limit(async () => {
      if (existing?.tree_sha === component.tree_sha && existing?.meta && existing?.git_sha) {
        component.meta = existing.meta;
        component.git_sha = existing.git_sha;
      } else {
        const [metaYml, commits] = await Promise.all([
          getGitHubFile("modules", component.path, commitSha),
          octokit.rest.repos.listCommits({
            owner: "nf-core",
            repo: "modules",
            path: component.path.replace("/meta.yml", ""),
            per_page: 1,
          }),
        ]);
        component.meta = metaYml ? parse(metaYml) : undefined;
        component.git_sha = commits.data[0]?.sha;
      }
      bar.tick();
    });

  // upsert items into target, replacing entries with the same name
  const mergeByName = (target, items) => {
    const indexByName = new Map(target.map((c, i) => [c.name, i]));
    for (const item of items) {
      const index = indexByName.get(item.name);
      if (index !== undefined) {
        target[index] = item;
      } else {
        target.push(item);
      }
    }
  };

  const modules = metaYmlFiles.filter((f) => f.path.startsWith("modules/")).map((f) => buildComponent(f, "module"));

  const existingModulesMap = new Map(components.modules.map((m) => [m.name, m]));
  let bar = new ProgressBar("  fetching module meta.ymls [:bar] :percent :etas", { total: modules.length });
  await Promise.all(modules.map((module) => enrichComponent(module, existingModulesMap.get(module.name), bar)));
  mergeByName(components.modules, modules);

  const subworkflows = metaYmlFiles
    .filter((f) => f.path.startsWith("subworkflows/") && !f.path.includes("homer/groseq"))
    .map((f) => buildComponent(f, "subworkflow"));

  const existingSubworkflowsMap = new Map((components.subworkflows ?? []).map((m) => [m.name, m]));
  bar = new ProgressBar("  fetching subworkflow meta.ymls [:bar] :percent :etas", { total: subworkflows.length });
  await Promise.all(
    subworkflows.map((subworkflow) => enrichComponent(subworkflow, existingSubworkflowsMap.get(subworkflow.name), bar)),
  );

  if (!components.subworkflows) {
    components.subworkflows = [];
  }
  mergeByName(components.subworkflows, subworkflows);

  // cross-reference: record on each (merged) module which subworkflows include it
  const modulesByName = new Map(components.modules.map((m) => [m.name, m]));
  for (const subworkflow of subworkflows) {
    for (const module of subworkflow.meta?.modules ?? []) {
      const existing = modulesByName.get(module);
      if (existing) {
        if (existing.subworkflows) {
          existing.subworkflows.push(subworkflow.name);
        } else {
          existing.subworkflows = [subworkflow.name];
        }
      }
    }
  }
  // Update pipelines that use modules and subworkflows
  let pipelinesToProcess = pipelines.remote_workflows;

  if (singlePipelineName) {
    pipelinesToProcess = pipelines.remote_workflows.filter((pipeline) => pipeline.name === singlePipelineName);
    if (pipelinesToProcess.length === 0) {
      console.error(`Pipeline ${singlePipelineName} not found in pipelines.json`);
      process.exit(1);
    }
    console.log(`Processing components for single pipeline: ${singlePipelineName}`);

    // Clean up existing pipeline references for this specific pipeline
    components.modules.forEach((module) => {
      if (module.pipelines) {
        module.pipelines = module.pipelines.filter((p) => p.name !== singlePipelineName);
        if (module.pipelines.length === 0) {
          delete module.pipelines;
        }
      }
    });

    components.subworkflows.forEach((subworkflow) => {
      if (subworkflow.pipelines) {
        subworkflow.pipelines = subworkflow.pipelines.filter((p) => p.name !== singlePipelineName);
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
              if (
                !components.subworkflows[index].pipelines.some(
                  (e) => e.name === entry.name && e.version === entry.version,
                )
              ) {
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
  components.commit_sha = commitSha;
  components.updated_at = new Date().toISOString();
  console.log("  writing components.json");
  writeFileSync(path.join(__dirname, "/public/components.json"), JSON.stringify(components, null, 2));
};

writeComponentsJson();
