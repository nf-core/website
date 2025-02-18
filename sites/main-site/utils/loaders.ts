// originally from https://github.com/gingerchew/astro-github-file-loader
import type { AstroConfig, MarkdownHeading } from "astro";
import type { Loader, LoaderContext } from "astro/loaders";
import { octokit } from "@components/octokit.js";

type GithubTreeLeaf = {
    path: string;
    mode: string;
    type: "tree" | "blob"; // tree is a directory, blob is a file
    sha: string;
    url: string;
};

type GithubTreeData = {
    url: string;
    hash: string;
    tree: GithubTreeLeaf[];
};
// Taken from astro content.d.ts
export interface RenderedContent {
    html: string;
    metadata?: {
        headings?: MarkdownHeading[];
        frontmatter?: Record<string, any>;
        imagePaths?: Array<string>;
        [key: string]: unknown;
    };
}

type ProcessorFileExtension = string;

type Processors = Record<ProcessorFileExtension, (str: string, config: AstroConfig) => Promise<RenderedContent>>;

interface PolicyLoaderConfig {
    org: string;
    repo: string;
    ref: string;
    path: string | RegExp; // Change path type to allow RegExp
    processors: Processors;
}

function createProcessors(processors: Processors) {
    return new Proxy(processors, {
        get(target, key: keyof Processors) {
            return key in target
                ? async (str: string, c: AstroConfig) => await target[key](str, c)
                : (_str: string, _c: AstroConfig) => ({
                      html: "",
                      metadata: {
                          error: `Could not find processor for extension: .${key}, are you sure you passed one in?`,
                      },
                  });
        },
    });
}

export function githubFileLoader({ org, repo, ref, processors, path }: PolicyLoaderConfig): Loader {
    const baseUrl = `https://raw.githubusercontent.com/${org}/${repo}/${ref}`;

    const get = async <T>(filepath: string, type: "json" | "text"): Promise<T> => {
        try {
            // If this is a GitHub API request, use octokit
            if (filepath.startsWith('https://api.github.com')) {
                const response = await octokit.request('GET ' + filepath.replace('https://api.github.com', ''));
                return response.data as T;
            }
            // Otherwise use fetch for raw content
            const response = await fetch(filepath);
            if (!response.ok) throw new Error(`HTTP ${response.status}`);
            return type === "json" ? await response.json() : await response.text();
        } catch (e) {
            console.error(`Failed to fetch ${filepath}:`, e);
            throw e;
        }
    };

    return {
        name: "github-file-loader",
        load: async ({ generateDigest, store, config, logger, meta }: LoaderContext) => {
            logger.info(`Loading files from ${org}/${repo}@${ref}`);

            // Get the last tree SHA we processed
            const lastTreeSha = meta.get('lastTreeSha');

            // Fetch current tree data using octokit directly
            const { data: treeData } = await octokit.rest.git.getTree({
                owner: org,
                repo: repo,
                tree_sha: ref,
                recursive: '1'
            });

            const currentTreeSha = treeData.sha;

            // If tree hasn't changed, we can skip processing
            if (lastTreeSha && lastTreeSha === currentTreeSha) {
                logger.info("No changes detected in GitHub tree, skipping update");
                return;
            }

            const $ = createProcessors(processors);

            // Track processed files to help clean up removed ones
            const processedIds = new Set<string>();

            for await (const leaf of treeData.tree) {
                // Skip directories and .github files
                if (leaf.type === "tree" || leaf.path?.includes(".github/")) continue;

                const pathMatch =
                    (path && typeof path === "string" && leaf.path?.includes(path)) ||
                    (path && typeof path !== "string" && path.test(leaf.path ?? ""));
                if (!pathMatch) continue;

                const [id, extension] = leaf.path?.split(".") ?? [];
                processedIds.add(id);

                // Check if we already have this file with the same SHA
                const existingEntry = store.get(id);
                if (existingEntry?.data?.sha === leaf.sha) {
                    logger.debug(`File ${leaf.path} unchanged, skipping`);
                    continue;
                }

                logger.info(`Processing ${leaf.path}`);

                const body = await get<string>(`${baseUrl}/${leaf.path}`, "text");
                const digest = generateDigest(body);

                const { html, metadata } = await $[extension as keyof Processors](body, config);

                // get the last commit date for the file using GitHub API
                const lastCommit = await get<any>(`https://api.github.com/repos/${org}/${repo}/commits?path=${leaf.path}&sha=${ref}`, "json");

                store.set({
                    id,
                    data: {
                        id,
                        extension,
                        org,
                        repo,
                        ref,
                        sha: leaf.sha, // Store SHA for future comparisons
                        lastCommit: lastCommit[0].commit.committer.date,
                    },
                    body,
                    rendered: { html, metadata },
                    digest
                });
            }

            // Remove entries that no longer exist in the tree
            for (const id of store.keys()) {
                if (!processedIds.has(id)) {
                    logger.info(`Removing deleted file: ${id}`);
                    store.delete(id);
                }
            }

            // Update the stored tree SHA
            meta.set('lastTreeSha', currentTreeSha);
            logger.info('GitHub file loader completed successfully');
        },
    };
}
