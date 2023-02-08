#! /usr/bin/env node
// const { promisify } = require('util');
// const { exec } = require('child_process');
// const { octokit } = require('../src/components/octokit');
// const readFile = promisify(fs.readFile);
// const writeFile = promisify(fs.writeFile);
import * as dotenv from 'dotenv';
import { readFileSync, writeFileSync } from 'fs';
import { Octokit } from 'octokit';
import path from 'path';

// see https://github.com/motdotla/dotenv#how-do-i-use-dotenv-with-import
dotenv.config();

const octokit = new Octokit({
  auth: process.env.GITHUB_TOKEN,
});

// get current path
const __dirname = path.resolve();

// write the pipelines.json file
const writePipelinesJson = async () => {
  const pipelinesJson = readFileSync(path.join(__dirname, '/public/pipelines.json'));
  const namesJson = readFileSync(path.join(__dirname, '/public/pipeline_names.json'));
  const pipelines = JSON.parse(pipelinesJson);
  const names = JSON.parse(namesJson).pipeline;
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
    // get the releases
    let { data: releases } = await octokit.rest.repos.listReleases({
      owner: 'nf-core',
      repo: name,
    });
    releases = [...releases, { tag_name: 'dev', published_at: Date.now() }];
    data['releases'] = releases.map(async (release) => {
      const { tag_name, published_at } = release;
      const doc_files = await octokit.rest.repos
        .getContent({
          owner: 'nf-core',
          repo: name,
          path: 'docs',
          ref: tag_name,
        })
        .then((response) => {
          return response.data
            .filter((file) => {
              return file.name.includes('.md');
            })
            .map((file) => {
              return file.path;
            });
        });
      return { tag_name, published_at, doc_files };
    });
    // resolve the promises
    data['releases'] = await Promise.all(data['releases']);
    // update in pipelines.remote_workflows if entry with name exists or add it otherwise
    const index = pipelines.remote_workflows.findIndex((workflow) => workflow.name === name);
    if (index > -1) {
      pipelines.remote_workflows[index] = data;
    } else {
      pipelines.remote_workflows.push(data);
    }

    // write the pipelines.json file
  }
  const json = JSON.stringify(pipelines, null, 4);
  await writeFileSync(path.join(__dirname, '/public/pipelines.json'), json, 'utf8');
};

writePipelinesJson();
