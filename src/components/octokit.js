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

export async function getCurrentRateLimitRemaining() {
  try {
    // Make a request to any endpoint (e.g., get the authenticated user)
    const response = await octokit.rest.repos.get({
      owner: 'nf-core',
      repo: 'nf-co.re',
    });

    // Get the 'x-ratelimit-remaining' header from the response headers
    const rateLimitRemaining = response.headers['x-ratelimit-remaining'];

    console.log(`Rate limit remaining: ${rateLimitRemaining}`);
  } catch (error) {
    console.error('Error occurred:', error);
  }
}

export const getGitHubFile = async (repo, path, ref) => {
  // console.log(`Getting ${path} from ${repo} ${ref}`);
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
        console.log(`something else happened for ${path} in ${repo} ${ref}`, error);
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
        // add github url to image links in markdown if they are relative
        content = content.replaceAll(/!\[(.*?)\]\((.*?)\)/g, (match, p1, p2) => {
          if (p2.startsWith('http')) {
            return match;
          } else {
            return `![${p1}](https://raw.githubusercontent.com/nf-core/${repo}/${ref}/${parent_directory}/${p2})`;
          }
        });
        // add github url to html img tags in markdown
        content = content.replaceAll(/<img(.*?)src="(.*?)"/g, (match, p1, p2) => {
          if (p2.startsWith('http')) {
            return match;
          } else {
            return `<img${p1}src="https://raw.githubusercontent.com/nf-core/${repo}/${ref}/${parent_directory}/${p2}"`;
          }
        });
        // add github url to html img tags in markdown for dark mode images
        content = content.replaceAll(/<source(.*?)srcset="(.*?)"/g, (match, p1, p2) => {
          if (p2.startsWith('http')) {
            return match;
          } else {
            return `<source${p1}src="https://raw.githubusercontent.com/nf-core/${repo}/${ref}/${parent_directory}/${p2}"`;
          }
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
  const getFilesInDir = async (directory) => {
    try {
      const response = await octokit.rest.repos.getContent({
        owner: 'nf-core',
        repo: pipeline,
        path: directory,
        ref: version,
      });
      const files = [];
      for (const file of response.data) {
        if (file.type === 'dir' && file.name !== 'images') {
          const subFiles = await getFilesInDir(file.path);
          files.push(...subFiles);
        } else if (
          file.type === 'file' &&
          file.name.includes('.md') &&
          (file.path.includes('output') || file.path.includes('usage'))
        ) {
          files.push(file.path);
        }
      }
      return files;
    } catch (error) {
      if (error.status === 404) {
        console.log(`Directory ${directory} not found in ${pipeline} ${version}`);
        console.log(error.request.url);
      } else {
        console.log(error);
      }
      return [];
    }
  };

  const doc_files = await getFilesInDir('docs');
  return doc_files;
};
