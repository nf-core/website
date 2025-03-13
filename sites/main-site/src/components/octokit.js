import * as dotenv from 'dotenv';
import { Octokit } from 'octokit';
import { retry } from "@octokit/plugin-retry";
import { throttling } from "@octokit/plugin-throttling";

// Get GitHub token from the appropriate environment
let GITHUB_TOKEN;

// Load environment variables
dotenv.config();

try {
  // Try to use Astro environment if available
  if (typeof import.meta !== 'undefined' && import.meta.env) {
    const astroEnv = await import('astro:env/server');
    GITHUB_TOKEN = astroEnv.GITHUB_TOKEN;
    console.log("Using GITHUB_TOKEN from Astro environment");
  } else {
    // Fall back to process.env
    GITHUB_TOKEN = process.env.GITHUB_TOKEN;
    console.log("Using GITHUB_TOKEN from process.env");
  }
} catch (error) {
  // Fall back to process.env if Astro import fails
  GITHUB_TOKEN = process.env.GITHUB_TOKEN;
  console.log("error:", error);
  console.log("Astro environment not detected, using process.env.GITHUB_TOKEN");
}

// For Astro, we need to access the token through import.meta.env.GITHUB_TOKEN
// Note that this needs to be available server-side
const MyOctokit = Octokit.plugin(retry, throttling);
export const octokit = new MyOctokit({
  auth: GITHUB_TOKEN,
  throttle: {
    onRateLimit: (retryAfter, options, octokit, retryCount) => {
      console.log(`Request quota exhausted for request ${options.method} ${options.url}. Retrying in ${retryAfter} seconds. Retry count: ${retryCount}`);

      if (retryCount < 1) {
        // only retries once
        console.log(`Retrying after ${retryAfter} seconds!`);
        return true;
      }
    },
    onSecondaryRateLimit: (retryAfter, options, octokit) => {
      // does not retry, only logs a warning
      console.log(`Secondary quota detected for request ${options.method} ${options.url}`);
    },
  },
  retry: {
    doNotRetry: ["429"],
  },
});

export default octokit;

export async function getCurrentRateLimitRemaining() {
  try {
    // Make a request to any endpoint (e.g., get the authenticated user)
    const { data: response } = await octokit.request('GET /rate_limit', {
      headers: {
        'X-GitHub-Api-Version': '2022-11-28'
      }
    });
    console.log(`Rate limit remaining: ${response.resources.core.remaining} out of ${response.resources.core.limit}, reset at ${new Date(response.resources.core.reset * 1000).toLocaleString()}`);
  } catch (error) {
    console.error('Error occurred:', error);
  }
}
export const githubFolderExists = async (repo, path, ref) => {
  return fetch(`https://github.com/nf-core/${repo}/tree/${ref}/${path}`).then(response  => response.ok);
}
export const getGitHubFile = async (repo, path, ref) => {
  try {
    const response = await fetch(`https://raw.githubusercontent.com/nf-core/${repo}/${ref}/${path}`);
    if (response.ok) {
      let content = await response.text();
      if (path.endsWith('.md') || path.endsWith('.mdx')) {
        const parent_directory = path.split('/').slice(0, -1).join('/');
        // add github url to image links in markdown if they are relative
        content = content.replaceAll(/!\[([^\][]*[?[^\][]*\]?[^\][]*)\]\((.*?)\)/g, (match, p1, p2) => {
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
        // prefix links to CONTRIBUTING.md, CITATIONS.md, CHANGELOG.md with github url
        content = content.replaceAll(
          /\[(.*?)\]\((\.github\/CONTRIBUTING\.md|CITATIONS\.md|CHANGELOG\.md)\)/g,
          (match, p1, p2) => {
            if (p2.startsWith('http')) {
              return match;
            } else {
              return `[${p1}](https://github.com/nf-core/${repo}/blob/${ref}/${p2})`;
            }
          },
        );
        // prefix links to files in the assets directory with github url
        content = content.replaceAll(/\[(.*?)\]\(((\.\.\/)*assets\/.*?)\)/g, (match, p1, p2) => {
          if (p2.startsWith('http')) {
            return match;
          } else {
            return `[${p1}](https://github.com/nf-core/${repo}/blob/${ref}/${p2.replace('../assets/', 'assets/')})`;
          }
        });

        // convert github style admonitions to docusaurus admonitions
        content = content.replace(
          /> \[!(NOTE|TIP|IMPORTANT|WARNING|CAUTION)\]\s*\n((?:> [^\n]*\s*?)+)/g,
          (match, type, content) => {
            const cleanedContent = content.replace(/> /g, '').trim();
            const admonitionType = type.toLowerCase();

            if (admonitionType === 'important') {
              return `:::info{title=Important}\n${cleanedContent}\n:::\n\n`;
            }
            if (admonitionType === 'caution') {
              return `:::danger{title=Caution}\n${cleanedContent}\n:::\n\n`;
            }

            return `:::${admonitionType}\n${cleanedContent}\n:::\n\n`;
          },
        );

        // remove .md(x) from links with anchor tags
        content = content.replaceAll(/\[([^\][]*)\]\((.*?)\.mdx?#(.*?)\)/g, '[$1]($2#$3)');

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
