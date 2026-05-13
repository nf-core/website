#! /usr/bin/env node
import { getGitHubFile, getCurrentRateLimitRemaining } from '../sites/main-site/src/components/octokit.js';
import { readFileSync, writeFileSync, existsSync } from 'fs';
import path from 'path';
import ProgressBar from 'progress';

const __dirname = path.resolve();

const args = process.argv.slice(2);
const singlePipelineName = args[0] || null;
const singleVersion = args[1] || null;

if (singleVersion && !singlePipelineName) {
  console.error('A pipeline name must be provided when specifying a version');
  process.exit(1);
}

console.log(await getCurrentRateLimitRemaining());

// Internal key combining param name + type so same-named params with different types stay separate.
const internalKey = (name, type) => `${name}|${type ?? ''}`;

// Convert the grouped output format { paramName: [{type, pipelines}] }
// back into the flat internal map { "paramName|type": {name, type, pipelines} }.
function ungroupParams(grouped) {
  const flat = {};
  for (const [name, variants] of Object.entries(grouped)) {
    for (const { type, pipelines } of variants) {
      flat[internalKey(name, type)] = { name, type, pipelines };
    }
  }
  return flat;
}

// Convert the flat internal map into the grouped output format.
function groupParams(flat) {
  const grouped = {};
  for (const { name, type, pipelines } of Object.values(flat)) {
    if (!grouped[name]) grouped[name] = [];
    grouped[name].push({ type, pipelines });
  }
  // Sort each variant list by type name, then sort params alphabetically by key.
  return Object.fromEntries(
    Object.keys(grouped)
      .sort()
      .map((name) => [
        name,
        grouped[name]
          .sort((a, b) => (a.type ?? '').localeCompare(b.type ?? ''))
          .map(({ type, pipelines }) => ({
            type,
            pipelines: pipelines.sort(
              (a, b) => a.name.localeCompare(b.name) || a.version.localeCompare(b.version),
            ),
          })),
      ]),
  );
}

export const writeParamsJson = async () => {
  const pipelinesJson = readFileSync(path.join(__dirname, '/public/pipelines.json'));
  const pipelines = JSON.parse(pipelinesJson);

  const paramsJsonPath = path.join(__dirname, '/public/params.json');

  // Internal map keyed by "paramName|type"
  let params = {};

  if (singlePipelineName) {
    if (existsSync(paramsJsonPath)) {
      params = ungroupParams(JSON.parse(readFileSync(paramsJsonPath)));
    }
    // Remove existing entries for the targeted pipeline+version before re-processing
    for (const key of Object.keys(params)) {
      params[key].pipelines = params[key].pipelines.filter((p) => {
        if (p.name !== singlePipelineName) return true;
        return singleVersion ? p.version !== singleVersion : false;
      });
      if (params[key].pipelines.length === 0) {
        delete params[key];
      }
    }

    const pipelineEntry = pipelines.remote_workflows.find((p) => p.name === singlePipelineName);
    if (!pipelineEntry) {
      console.error(`Pipeline ${singlePipelineName} not found in pipelines.json`);
      process.exit(1);
    }
    if (singleVersion && !pipelineEntry.releases.some((r) => r.tag_name === singleVersion)) {
      console.error(`Version ${singleVersion} not found for pipeline ${singlePipelineName} in pipelines.json`);
      process.exit(1);
    }
    console.log(`Processing parameters for: ${singlePipelineName}${singleVersion ? `@${singleVersion}` : ''}`);
  }

  const pipelinesToProcess = singlePipelineName
    ? pipelines.remote_workflows.filter((p) => p.name === singlePipelineName)
    : pipelines.remote_workflows;

  // Build a set of already-processed {pipeline@version} pairs so re-runs skip unchanged releases
  const processed = new Set();
  if (!singlePipelineName) {
    for (const { pipelines: pipelineList } of Object.values(params)) {
      for (const { name, version } of pipelineList) {
        processed.add(`${name}@${version}`);
      }
    }
  }

  const CONCURRENCY_LIMIT = 10;
  const bar = new ProgressBar('  fetching schemas [:bar] :percent :etas', { total: pipelinesToProcess.length });

  const processPipeline = async (pipeline) => {
    const targetRelease = singleVersion
      ? pipeline.releases?.find((r) => r.tag_name === singleVersion && r.has_schema)
      : (pipeline.releases?.find((r) => r.tag_name !== 'dev' && r.has_schema) ??
         pipeline.releases?.find((r) => r.tag_name === 'dev' && r.has_schema));

    if (!targetRelease) {
      if (singleVersion) {
        console.error(`Version ${singleVersion} not found or has no schema for ${pipeline.name}`);
      }
      bar.tick();
      return;
    }

    // Skip if this exact release was already parsed into params.json
    if (processed.has(`${pipeline.name}@${targetRelease.tag_name}`)) {
      bar.tick();
      return;
    }

    const schemaContent = await getGitHubFile(pipeline.name, 'nextflow_schema.json', targetRelease.tag_name);

    if (!schemaContent) {
      bar.tick();
      return;
    }

    let schema;
    try {
      schema = JSON.parse(schemaContent);
    } catch {
      bar.tick();
      return;
    }

    // Collect param names + types from definitions (draft-07) or $defs (draft 2020-12)
    const defs = schema.definitions || schema['$defs'] || {};
    const paramMap = new Map(); // paramName -> type

    for (const def of Object.values(defs)) {
      if (def.properties) {
        for (const [paramName, propValue] of Object.entries(def.properties)) {
          const t = propValue.type ?? null;
          paramMap.set(paramName, Array.isArray(t) ? t.filter((x) => x !== 'null')[0] ?? t[0] ?? null : t);
        }
      }
    }

    // Top-level properties that aren't $ref pointers to definition groups
    if (schema.properties) {
      for (const [paramName, propValue] of Object.entries(schema.properties)) {
        if (!propValue.$ref) {
          const t = propValue.type ?? null;
          paramMap.set(paramName, Array.isArray(t) ? t.filter((x) => x !== 'null')[0] ?? t[0] ?? null : t);
        }
      }
    }

    const entry = { name: pipeline.name, version: targetRelease.tag_name };
    for (const [paramName, type] of paramMap) {
      const key = internalKey(paramName, type);
      if (!params[key]) {
        params[key] = { name: paramName, type, pipelines: [] };
      }
      params[key].pipelines.push(entry);
    }

    bar.tick();
  };

  for (let i = 0; i < pipelinesToProcess.length; i += CONCURRENCY_LIMIT) {
    const batch = pipelinesToProcess.slice(i, i + CONCURRENCY_LIMIT);
    await Promise.all(batch.map(processPipeline));
  }

  console.log('  writing params.json');
  writeFileSync(paramsJsonPath, JSON.stringify(groupParams(params), null, 2));
};

writeParamsJson();
