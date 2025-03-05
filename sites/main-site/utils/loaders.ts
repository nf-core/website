// originally from https://github.com/gingerchew/astro-github-file-loader
import type { AstroConfig, MarkdownHeading } from "astro";
import type { Loader, LoaderContext } from "astro/loaders";
import { octokit, getCurrentRateLimitRemaining } from "@components/octokit.js";
import { z } from "zod";

import { createMarkdownProcessor } from "@astrojs/markdown-remark";
import type { MarkdownProcessor } from "@astrojs/markdown-remark";
// Import the remark plugin directly
import remarkGitHubMarkdown from "../../../bin/remark-github-markdown.js";

type GithubTreeLeaf = {
    path: string;
    mode: string;
    type: "tree" | "blob"; // tree is a directory, blob is a file
    sha: string;
    url: string;
};

type GitHubCommit = {
    commit: {
        committer: {
            date: string;
        };
    };
}[];

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

// Update the type definition to include the optional filepath parameter
type Processors = Record<
    ProcessorFileExtension,
    (str: string, config: AstroConfig, filepath?: string) => Promise<RenderedContent>
>;

interface PolicyLoaderConfig {
    org: string;
    repo: string;
    ref: string;
    path: string | RegExp; // Change path type to allow RegExp
    processors: Processors;
    getDates: boolean;
}

function createProcessors(processors: Processors) {
    return new Proxy(processors, {
        get(target, key: keyof Processors) {
            return key in target
                ? async (str: string, c: AstroConfig, filepath?: string) => await target[key](str, c, filepath)
                : (_str: string, _c: AstroConfig, _filepath?: string) => ({
                      html: "",
                      metadata: {
                          error: `Could not find processor for extension: .${key}, are you sure you passed one in?`,
                      },
                  });
        },
    });
}

const mdProcessors = new Map<AstroConfig, MarkdownProcessor>();

// Simplified markdown processor that delegates all transformations to the remark plugin
function createGitHubMarkdownProcessor(repo: string, ref: string) {
    return {
        md: async (text: string, config: AstroConfig, filepath?: string): Promise<RenderedContent> => {
            // Extract parent directory from the filepath if provided
            const parent_directory = filepath
                ? filepath.includes("/")
                    ? filepath.split("/").slice(0, -1).join("/")
                    : ""
                : "";

            // Get or create the markdown processor
            const processor = (
                mdProcessors.has(config)
                    ? mdProcessors.get(config)
                    : mdProcessors
                          .set(
                              config,
                              await createMarkdownProcessor({
                                  ...config.markdown,
                                  remarkPlugins: [
                                      [
                                          remarkGitHubMarkdown,
                                          {
                                              repo,
                                              ref,
                                              parent_directory,
                                          },
                                      ],
                                      ...(config.markdown.remarkPlugins || []),
                                      // Use the imported function reference instead of a string
                                  ],
                              }),
                          )
                          .get(config)
            )!;

            try {
                // Render the markdown
                const { code: html, metadata } = await processor.render(text);

                return { html, metadata };
            } catch (error) {
                console.error(`Error processing markdown for ${repo}/${ref}:`, error);

                // Provide a fallback with error information
                return {
                    html: `<div class="error">Error processing markdown: ${error.message}</div>`,
                    metadata: {},
                };
            }
        },
    };
}

class GitHubContentFetcher {
    private baseUrl: string;

    constructor(
        private org: string,
        private repo: string,
        private ref: string,
        private context: LoaderContext,
    ) {
        this.baseUrl = `https://raw.githubusercontent.com/${org}/${repo}/${ref}`;
    }

    private async get<T>(filepath: string, type: "json" | "text"): Promise<T> {
        try {
            if (filepath.startsWith("https://api.github.com")) {
                const response = await octokit.request("GET " + filepath.replace("https://api.github.com", ""));
                return response.data as T;
            }
            const response = await fetch(filepath);
            if (!response.ok) throw new Error(`HTTP ${response.status}`);
            return type === "json" ? await response.json() : await response.text();
        } catch (e) {
            console.error(`Failed to fetch ${filepath}:`, e);
            throw e;
        }
    }

    async processFile(
        filepath: string,
        processors: Processors,
        getDates: boolean,
        additionalData: Record<string, unknown> = {},
    ) {
        const { store, generateDigest, config, logger } = this.context;
        logger.debug(`Processing ${this.org}/${this.repo}@${this.ref}/${filepath}`);

        let [id, extension] = filepath.split(".") ?? [];
        id = this.repo + "/" + this.ref + "/" + id;

        // check if the file exists in the store
        const existingEntry = store.get(id);
        if (existingEntry && existingEntry.data.repo === this.repo && existingEntry.data.ref === this.ref) {
            logger.debug(`File ${filepath} unchanged, skipping`);
            return;
        }

        let body = await this.get<string>(`${this.baseUrl}/${filepath}`, "text");
        const digest = generateDigest(body);

        // Pass repo and ref information to the markdown processor via config
        const processorConfig = {
            ...config,
            data: {
                repo: this.repo,
                ref: this.ref,
                parent_directory: filepath.includes("/") ? filepath.split("/").slice(0, -1).join("/") : "",
            },
        };

        try {
            logger.debug(`Processing ${filepath}, content length: ${body.length}`);
            // Pass the filepath as the third argument
            const { html, metadata } = await createProcessors(processors)[extension as keyof Processors](
                body,
                processorConfig,
                filepath,
            );
            logger.debug(`Successfully processed ${filepath}`);

            let lastCommit: GitHubCommit | null = null;
            if (getDates) {
                lastCommit = await this.get<GitHubCommit>(
                    `https://api.github.com/repos/${this.org}/${this.repo}/commits?path=${filepath}&ref=${this.ref}`,
                    "json",
                );
            }

            logger.debug(`Setting ${id} in store`);
            store.set({
                id,
                data: {
                    id,
                    extension,
                    org: this.org,
                    repo: this.repo,
                    ref: this.ref,
                    lastCommit: lastCommit?.[0].commit.committer.date || null,
                    ...additionalData,
                },
                body,
                rendered: { html, metadata },
                digest,
            });

            return id;
        } catch (error) {
            logger.error(`Error processing file ${filepath}: ${error.message}`);
            logger.error(error.stack);
        }
    }

    async processFilesFromTree({
        path,
        processors,
        getDates = false,
    }: {
        path: string | RegExp;
        processors: Processors;
        getDates?: boolean;
    }) {
        const { logger, meta, store } = this.context;
        logger.debug(`Loading files from ${this.org}/${this.repo}@${this.ref}`);

        const lastTreeSha = meta.get("lastTreeSha");
        const treeData = await this.get<{ sha: string; tree: GithubTreeLeaf[] }>(
            `https://api.github.com/repos/${this.org}/${this.repo}/git/trees/${this.ref}?recursive=1`,
            "json",
        );
        const currentTreeSha = treeData.sha;

        if (lastTreeSha && lastTreeSha === currentTreeSha) {
            logger.debug("No changes detected in GitHub tree, skipping update");
            return;
        }

        const processedIds = new Set<string>();

        for await (const leaf of treeData.tree) {
            if (leaf.type === "tree" || leaf.path?.includes(".github/")) continue;

            const pathMatch =
                (path && typeof path === "string" && leaf.path?.includes(path)) ||
                (path && typeof path !== "string" && path.test(leaf.path ?? ""));
            if (!pathMatch) continue;

            logger.debug(`Processing ${leaf.path}`);
            const id = await this.processFile(leaf.path, processors, getDates);
            if (id) {
                processedIds.add(id);
            }
        }

        // Cleanup and update tree SHA
        for (const id of store.keys()) {
            if (!processedIds.has(id)) {
                logger.debug(`Removing deleted file: ${id}`);
                store.delete(id);
            }
        }

        meta.set("lastTreeSha", currentTreeSha);
        logger.debug("GitHub file processing completed successfully");
    }
}

export function githubFileLoader({ org, repo, ref, processors, path, getDates = false }: PolicyLoaderConfig): Loader {
    return {
        name: "github-file-loader",
        load: async (context: LoaderContext) => {
            const fetcher = new GitHubContentFetcher(org, repo, ref, context);
            await fetcher.processFilesFromTree({ path, processors, getDates });
        },
    };
}

const schemaProcessor = {
    json: async (text: string, config: AstroConfig): Promise<RenderedContent> => {
        const schema = JSON.parse(text);
        const definitions = schema.definitions || schema.$defs || schema.properties || {};
        let headings: { slug: string; text: string; depth: number; fa_icon: string; hidden: boolean }[] = [];
        if (Object.keys(definitions).length > 0) {
            headings = Object.entries(definitions).map(([key, value]: [string, any]) => ({
                slug: key.replaceAll("_", "-"),
                text: value?.title || key,
                depth: 1,
                fa_icon: value?.fa_icon || "",
                hidden: value?.properties && Object.values(value.properties).every((prop: any) => prop?.hidden),
            }));
        }
        return {
            html: text,
            metadata: {
                schema: schema["$schema"],
                headings,
            },
        };
    },
};

export function pipelineLoader(pipelines_json: {
    remote_workflows: {
        name: string;
        releases: { tag_name: string; doc_files: string[]; has_schema: boolean }[];
        archived: boolean;
        description: string;
        topics: string[];
        stargazers_count: number;
    }[];
}): Loader {
    return {
        name: "pipeline-loader",
        load: async (context: LoaderContext) => {
            console.log('Starting pipeline loader');
            getCurrentRateLimitRemaining();
            // Process pipelines sequentially
            for (const pipeline of pipelines_json.remote_workflows) {
                // Process releases sequentially
                for (const release of pipeline.releases) {
                    const fetcher = new GitHubContentFetcher("nf-core", pipeline.name, release.tag_name, context);
                    if (release.has_schema) {
                        // load the schema
                        await fetcher.processFile(`nextflow_schema.json`, schemaProcessor, false, {
                            name: pipeline.name,
                        });
                    }
                    release.doc_files.push("README.md");

                    // Create processor and fetcher once per pipeline/release
                    const processors = createGitHubMarkdownProcessor(pipeline.name, release.tag_name);

                    // Common metadata for all files in this pipeline/release
                    const metadata = {
                        name: pipeline.name,
                        archived: pipeline.archived,
                        releases: pipeline.releases,
                        description: pipeline.description,
                        topics: pipeline.topics,
                        stargazers_count: pipeline.stargazers_count,
                    };

                    // Process files sequentially
                    for (const doc_file of release.doc_files) {
                        context.logger.debug(`Loading ${pipeline.name}@${release.tag_name}/${doc_file}`);
                        try {
                            await fetcher.processFile(doc_file, processors, false, metadata);
                        } catch (error) {
                            context.logger.error(
                                `Error processing ${pipeline.name}@${release.tag_name}/${doc_file}: ${error.message}`,
                            );
                            // Continue with next file instead of failing the entire pipeline
                        }
                    }
                }
            }
            console.log('Pipeline loader completed successfully');
        },
    };
}

export function releaseLoader(pipelines_json: {
    remote_workflows: {
        name: string;
        releases: { tag_name: string; doc_files: string[]; has_schema: boolean }[];
        archived: boolean;
        description: string;
        topics: string[];
        stargazers_count: number;
    }[];
}): Loader {
    return {
        name: "release-notes-loader",
        load: async (context: LoaderContext) => {
            console.log('Starting release notes loader');
            getCurrentRateLimitRemaining();
            const { store, generateDigest, config, logger } = context;

            // Process pipelines sequentially
            for (const pipeline of pipelines_json.remote_workflows) {
                logger.debug(`Loading release notes for ${pipeline.name}`);

                try {
                    // Fetch all releases for this pipeline
                    const releases = await octokit.request(`GET /repos/nf-core/${pipeline.name}/releases`);

                    // Process each release
                    for (const release of releases.data) {
                        const releaseId = `${pipeline.name}/${release.tag_name}/release-notes`;

                        // Check if we already have this release in the store
                        const existingEntry = store.get(releaseId);
                        if (existingEntry && existingEntry.data.tag_name === release.tag_name) {
                            logger.debug(`Release notes for ${pipeline.name}@${release.tag_name} unchanged, skipping`);
                            continue;
                        }

                        // Get the release body (markdown content)
                        const body = release.body || '';
                        const digest = generateDigest(body);

                        // Create a markdown processor for this content
                        const processor = await createMarkdownProcessor(config.markdown);

                        try {
                            // Render the markdown
                            const { code: html, metadata } = await processor.render(body);

                            // Store the processed release notes
                            store.set({
                                id: releaseId,
                                data: {
                                    id: releaseId,
                                    name: pipeline.name,
                                    tag_name: release.tag_name,
                                    published_at: release.published_at,
                                    html_url: release.html_url,
                                    prerelease: release.prerelease,
                                    draft: release.draft,
                                    author: release.author ? {
                                        login: release.author.login,
                                        avatar_url: release.author.avatar_url,
                                        html_url: release.author.html_url
                                    } : null,
                                    // Include additional pipeline metadata
                                    archived: pipeline.archived,
                                    description: pipeline.description,
                                    topics: pipeline.topics,
                                    stargazers_count: pipeline.stargazers_count
                                },
                                body,
                                rendered: { html, metadata },
                                digest
                            });

                            logger.debug(`Successfully processed release notes for ${pipeline.name}@${release.tag_name}`);
                        } catch (error) {
                            logger.error(`Error processing release notes for ${pipeline.name}@${release.tag_name}: ${error.message}`);
                            // Continue with next release instead of failing the entire pipeline
                        }
                    }
                } catch (error) {
                    logger.error(`Error fetching releases for ${pipeline.name}: ${error.message}`);
                    // Continue with next pipeline instead of failing the entire process
                }

            }

            logger.info("GitHub release notes processing completed successfully");
        }
    };
}

// Export a default markdown processor for use in other files
export const md = {
    md: async (text: string, config: AstroConfig, filepath?: string): Promise<RenderedContent> => {
        const processor = await createMarkdownProcessor(config.markdown);
        try {
            const { code: html, metadata } = await processor.render(text);
            return { html, metadata };
        } catch (error) {
            console.error(`Error processing markdown: ${error.message}`);
            throw error;
        }
    },
};
