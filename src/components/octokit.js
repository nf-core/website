import * as dotenv from 'dotenv';
import { Octokit } from 'octokit';


if (!import.meta.env) {
  dotenv.config();
}
const octokit = new Octokit({
  // different env vars for node (used for scripts in `bin/`) and astro
  auth: import.meta.env ? import.meta.env.GITHUB_TOKEN : process.env.GITHUB_TOKEN,
});
export default octokit;

export const getGitHubFile = async (repo, path, ref) => {
  const response = await octokit
    .request('GET /repos/nf-core/{repo}/contents/{path}?ref={ref}', {
      repo: repo,
      path: path,
      ref: ref,
    })
    .catch((error) => {
      if (error.status === 404) {
        console.log(`File ${path} not found in ${repo} ${ref}`);
        console.log(error.request.url);
        return;
      } else {
        console.log(error);
        return;
      }
    })
    .then((response) => {
      if (response == null) {
        return;
      }
      let content = Buffer.from(response.data.content, 'base64').toString('utf-8');
      if (path.endsWith('.md')) {
        const parent_directory = path.split('/').slice(0, -1).join('/');
        // add github url to image links in markdown
        content = content.replaceAll(/!\[(.*?)\]\((.*?)\)/g, (match, p1, p2) => {
          return `![${p1}](https://raw.githubusercontent.com/nf-core/${repo}/${ref}/${parent_directory}/${p2})`;
        });
        // add github url to html img tags in markdown
        content = content.replaceAll(/<img(.*?)src="(.*?)"/g, (match, p1, p2) => {
          return `<img${p1}src="https://raw.githubusercontent.com/nf-core/${repo}/${ref}/${parent_directory}/${p2}"`;
        });
        // remove github warning and everything before from docs
        content = content.replace(/(.*?)(## :warning:)(.*?)(f)/s, '');
        // remove blockquote ending in "files._" from the start of the document
        content = content.replace(/(.*?)(files\._)/s, '');
        // cleanup heading
        content = content.replace('# nf-core/' + repo + ': ', '# ');
        // remove everything before introduction
        content = content.replace(/.*?## Introduction/gs, '## Introduction');
      }
      return content;
    });
  return response;
};

export const getDocFiles = async (pipeline, version) => {
  const doc_files = await octokit.rest.repos
    .getContent({
      owner: 'nf-core',
      repo: pipeline,
      path: 'docs',
      ref: version,
    })
    .then((response) => {
      return response.data
        .filter((file) => {
          return file.name.includes('.md') && !file.name.includes('README');
        })
        .map((file) => {
          return file.path;
        });
    });
  return doc_files;
};