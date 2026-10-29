// Walltime benchmarks for full site builds, reported to CodSpeed.
// Usage: node bench/build.bench.mjs [site ...]   (defaults to main-site and pipelines)
// Outside the CodSpeed runner this falls back to plain tinybench and prints a table.
import { Bench } from "tinybench";
import { withCodSpeed } from "@codspeed/tinybench-plugin";
import { spawnSync } from "node:child_process";
import { readFileSync, rmSync, writeFileSync } from "node:fs";
import path from "node:path";

const root = path.resolve(import.meta.dirname, "..");
const args = process.argv.slice(2);
const sites = args.length ? args : ["main-site", "pipelines"];

const pipelinesJsonPath = path.join(root, "public", "pipelines.json");
const pipelinesJson = readFileSync(pipelinesJsonPath, "utf8");

// Mimic the automated "Update json files" commit that precedes most deploys, where a few repo stats change
function touchPipelinesJson() {
    const data = JSON.parse(pipelinesJson);
    const pipeline = data.remote_workflows[0];
    pipeline.stargazers_count += 1;
    pipeline.watchers_count += 1;
    pipeline.watchers += 1;
    writeFileSync(pipelinesJsonPath, JSON.stringify(data, null, 4));
}

function restorePipelinesJson() {
    writeFileSync(pipelinesJsonPath, pipelinesJson);
}

function build(site) {
    const result = spawnSync("npm", ["run", "build", "-w", `sites/${site}`], { cwd: root, stdio: "inherit" });
    if (result.status !== 0) {
        throw new Error(`build of sites/${site} failed with exit code ${result.status}`);
    }
}

function clean(site) {
    // cacheDir as set in astro.config.base.mjs: incremental build manifest, content layer data store and fonts
    rmSync(path.join(root, "node_modules", ".astro", site), { recursive: true, force: true });
    rmSync(path.join(root, "sites", site, "dist"), { recursive: true, force: true });
}

// Builds take minutes, so run each task a fixed number of times instead of filling a time budget
const bench = withCodSpeed(new Bench({ iterations: 1, time: 0, warmup: false, throws: true }));
// async: false stops tinybench from calling the task once up front to detect a promise (an extra full build)
const taskOptions = { async: false };

for (const site of sites) {
    bench.add(`build ${site} (cold)`, () => build(site), { ...taskOptions, beforeEach: () => clean(site) });
    // Tasks run in order, so the cold task leaves a filled cache behind. Changing pipelines.json on top of it
    // measures experimental.incrementalBuild the way a typical deploy hits it, not the no-change best case
    bench.add(`build ${site} (warm)`, () => build(site), {
        ...taskOptions,
        beforeEach: touchPipelinesJson,
        afterEach: restorePipelinesJson,
    });
}

await bench.run();
console.table(bench.table());
