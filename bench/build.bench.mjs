// Walltime benchmarks for full site builds, reported to CodSpeed.
// Usage: node bench/build.bench.mjs [site ...]   (defaults to main-site and pipelines)
// Outside the CodSpeed runner this falls back to plain tinybench and prints a table.
import { Bench } from "tinybench";
import { withCodSpeed } from "@codspeed/tinybench-plugin";
import { spawnSync } from "node:child_process";
import { rmSync } from "node:fs";
import path from "node:path";

const root = path.resolve(import.meta.dirname, "..");
const args = process.argv.slice(2);
const sites = args.length ? args : ["main-site", "pipelines"];

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
const bench = withCodSpeed(new Bench({ iterations: 2, time: 0, warmup: false, throws: true }));
// async: false stops tinybench from calling the task once up front to detect a promise (an extra full build)
const taskOptions = { async: false };

for (const site of sites) {
    bench.add(`build ${site} (cold)`, () => build(site), { ...taskOptions, beforeEach: () => clean(site) });
    // Tasks run in order, so the cold task leaves a filled cache behind. With no source changes between
    // runs, this measures the best case for experimental.incrementalBuild
    bench.add(`build ${site} (warm)`, () => build(site), taskOptions);
}

await bench.run();
console.table(bench.table());
