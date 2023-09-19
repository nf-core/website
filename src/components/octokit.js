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
  try {
    const response = await fetch(`https://raw.githubusercontent.com/nf-core/${repo}/${ref}/${path}`);
    if (response.ok) {
      let content = await response.text();
      if (path.endsWith('.md') || path.endsWith('.mdx')) {
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
        // prefix links to CONTRIBUTING.md and CITATIONS.md with github url
        content = content.replaceAll(/\[(.*?)\]\((\.github\/CONTRIBUTING\.md|CITATIONS\.md)\)/g, (match, p1, p2) => {
          if (p2.startsWith('http')) {
            return match;
          } else {
            return `[${p1}](https://github.com/nf-core/${repo}/blob/${ref}/${p2})`;
          }
        });
        // remove github warning and everything before from docs
        content = content.replace(/(.*?)(## :warning:)(.*?)usage\)/s, '');
        // remove blockquote ending in "files._" from the start of the document
        content = content.replace(/(.*?)(files\._)/s, '');
        // cleanup heading
        content = content.replace('# nf-core/' + repo + ': ', '# ');
        // remove everything before introduction
        content = content.replace(/.*?## Introduction/gs, '## Introduction');
        // replace nextflow with groovy code blocks TODO: remove this when we have a nextflow syntax highlighter works in shiki
        content = content.replace(/```nextflow/g, '```groovy');
      }
      return content;
    } else {
      // console.log(`File ${path} not found in ${repo} ${ref}`);
      console.log(response.url);
      return null;
    }
  } catch (error) {
    if (error.status === 404) {
      console.log(`File ${path} not found in ${repo} ${ref}`);
      console.log(error.request.url);
    } else {
      console.log(error);
    }
    return null;
  }
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
