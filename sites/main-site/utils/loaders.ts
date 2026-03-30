// originally from https://github.com/gingerchew/astro-github-file-loader
import type { AstroConfig, MarkdownHeading } from "astro";
import type { Loader, LoaderContext } from "astro/loaders";
import { octokit, getCurrentRateLimitRemaining } from "@components/octokit.js";
import ProgressBar from "progress";
import { S3Client, ListObjectsV2Command } from "@aws-sdk/client-s3";
import { createMarkdownProcessor } from "@astrojs/markdown-remark";
import remarkGitHubMarkdown from "../../../bin/remark-github-markdown.js";

// ========================================
// TYPE DEFINITIONS
// ========================================

interface PipelineJson {
    remote_workflows: PipelineWorkflow[];
}

interface PipelineWorkflow {
    name: string;
    releases: PipelineRelease[];
    archived: boolean;
    description: string;
    topics: string[];
    stargazers_count: number;
}

interface PipelineRelease {
    tag_name: string;
    doc_files: string[];
    has_schema: boolean;
    tag_sha: string;
}

interface GithubTreeLeaf {
    path: string;
    mode: string;
    type: "tree" | "blob";
    sha: string;
    url: string;
}

interface GitHubCommit {
    commit: {
        committer: {
            date: string;
        };
    };
}

interface RenderedContent {
    html: string;
    metadata?: {
        headings?: MarkdownHeading[];
        frontmatter?: Record<string, any>;
        imagePaths?: Array<string>;
        [key: string]: unknown;
    };
}

type ProcessorFunction = (str: string, config: AstroConfig, filepath?: string) => Promise<RenderedContent>;
type Processors = Record<string, ProcessorFunction>;

interface GithubLoaderConfig {
    org: string;
    repo: string;
    ref: string;
    path: string | RegExp;
    processors: Processors;
    getDates?: boolean;
}

interface S3Object {
    Key: string;
    LastModified: Date;
    ETag: string;
    Size: number;
    StorageClass: string;
}

interface S3Prefix {
    Prefix: string;
}

// ========================================
// CONCURRENCY HELPERS
// ========================================

/**
 * Runs async tasks with a maximum concurrency limit
 */
async function parallelLimit<T>(tasks: (() => Promise<T>)[], limit: number): Promise<T[]> {
    const results: (T | undefined)[] = new Array(tasks.length);
    const queue = tasks.map((task, i) => ({ task, i }));
    let next = 0;
    const workers = Array.from({ length: Math.min(limit, tasks.length) }, async () => {
        while (next < queue.length) {
            const { task, i } = queue[next++];
            results[i] = await task();
        }
    });
    await Promise.all(workers);
    return results as T[];
}

// ========================================
// CONTENT PROCESSORS
// ========================================

/**
 * Creates a proxy for processors with fallback handling
 */
function createProcessors(processors: Processors) {
    return new Proxy(processors, {
        get(target, key: string) {
            return key in target
                ? async (str: string, c: AstroConfig, filepath?: string) => await target[key](str, c, filepath)
                : (_str: string, _c: AstroConfig, _filepath?: string) => ({
                      html: "",
                      metadata: {
                          error: `Could not find processor for extension: .${key}. Are you sure you passed one in?`,
                      },
                  });
        },
    });
}

/**
 * Creates a markdown processor with GitHub-specific transformations
 */
function createGitHubMarkdownProcessor(renderMarkdown: Function, options: { repo: string; ref: string }): Processors {
    return {
        md: async (text: string, config: AstroConfig, filepath?: string): Promise<RenderedContent> => {
            try {
                const parentDir = filepath?.includes("/") ? filepath.split("/").slice(0, -1).join("/") : "";
                const processor = await createMarkdownProcessor({
                    ...config.markdown,
                    remarkPlugins: [
                        [
                            remarkGitHubMarkdown,
                            {
                                org: "nf-core",
                                repo: options.repo,
                                ref: options.ref,
                                parent_directory: parentDir,
                            },
                        ],
                        ...(config.markdown?.remarkPlugins || []),
                    ],
                });
                const { code: html, metadata } = await processor.render(text);
                return { html, metadata };
            } catch (error) {
                console.error(`Error processing markdown for ${options.repo}/${options.ref}:`, error);
                return {
                    html: `<div class="error">Error processing markdown: ${error.message}</div>`,
                    metadata: {},
                };
            }
        },
    };
}

/**
 * Creates a schema processor for JSON schema files
 */
function createSchemaProcessor(renderMarkdown: Function): Processors {
    return {
        json: async (text: string, config: AstroConfig, filepath?: string): Promise<RenderedContent> => {
            try {
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

                    // Process markdown content in schema definitions
                    const defsKey = schema.definitions ? "definitions" : "$defs";
                    if (schema[defsKey]) {
                        // Deep copy the definitions to avoid mutation
                        schema[defsKey] = JSON.parse(JSON.stringify(schema[defsKey]));

                        await Promise.all(
                            Object.entries(schema[defsKey]).map(async ([key, definition]) => {
                                const def = definition as {
                                    description: string;
                                    help_text: string;
                                    properties: Record<string, any>;
                                };

                                // Render markdown for definition-level content
                                await Promise.all([
                                    renderMarkdownField(def, "description", renderMarkdown, `${key}`),
                                    renderMarkdownField(def, "help_text", renderMarkdown, `${key}`),
                                ]);

                                // Render markdown for property-level content
                                if (def.properties) {
                                    await Promise.all(
                                        Object.entries(def.properties).map(async ([propKey, property]) => {
                                            const prop = property as any;
                                            await Promise.all([
                                                renderMarkdownField(
                                                    prop,
                                                    "description",
                                                    renderMarkdown,
                                                    `${key}.${propKey}`,
                                                ),
                                                renderMarkdownField(
                                                    prop,
                                                    "help_text",
                                                    renderMarkdown,
                                                    `${key}.${propKey}`,
                                                ),
                                            ]);
                                        }),
                                    );
                                }
                            }),
                        );
                    }
                }

                return {
                    html: JSON.stringify(schema, null, 2),
                    metadata: {
                        schema: schema["$schema"],
                        headings,
                    },
                };
            } catch (error) {
                console.error(`Error processing schema for ${filepath}:`, error.message);
                return {
                    html: `<div class="error">Error processing schema: ${error.message}</div>`,
                    metadata: { error: error.message },
                };
            }
        },
    };
}

/**
 * Helper function to render markdown fields in schema objects
 */
async function renderMarkdownField(
    obj: any,
    fieldName: string,
    renderMarkdown: Function,
    context: string,
): Promise<void> {
    if (obj[fieldName]) {
        try {
            const { html } = await renderMarkdown(obj[fieldName]);
            obj[`${fieldName}_rendered`] = html;
        } catch (error) {
            console.warn(`Failed to render ${fieldName} for ${context}:`, error.message);
        }
    }
}

// ========================================
// GITHUB CONTENT FETCHER
// ========================================

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

    /**
     * Generic method to fetch data from GitHub
     */
    private async fetchData<T extends any>(filepath: string, type: "json" | "text", retries = 3): Promise<T> {
        try {
            if (filepath.startsWith("https://api.github.com")) {
                const response = await octokit.request("GET " + filepath.replace("https://api.github.com", ""));
                return response.data as T;
            }

            const response = await fetch(filepath);
            if (response.status === 429 && retries > 0) {
                const retryAfter = Number(response.headers.get("retry-after") ?? 60);
                console.log(`[rate-limit] ${filepath} — retrying in ${retryAfter}s (${retries} attempts left)`);
                await new Promise((resolve) => setTimeout(resolve, retryAfter * 1000));
                return this.fetchData<T>(filepath, type, retries - 1);
            }
            if (!response.ok) throw new Error(`HTTP ${response.status}`);

            return (type === "json" ? await response.json() : await response.text()) as T;
        } catch (error) {
            console.error(`Failed to fetch ${filepath}:`, error);
            throw error;
        }
    }

    /**
     * Process a single file from the repository
     */
    async processFile(
        filepath: string,
        processors: Processors,
        getDates: boolean,
        additionalData: Record<string, unknown> = {},
    ): Promise<string | undefined> {
        const { store, generateDigest, config, logger, renderMarkdown } = this.context;
        logger.info(`Processing ${this.org}/${this.repo}@${this.ref}/${filepath}`);

        // Generate unique ID for this file
        const [id, extension] = filepath.split(".") ?? [];
        const fileId = `${this.repo}/${this.ref}/${id}`;

        // Skip if file hasn't changed (except for configs)
        const existingEntry = store.get(fileId);
        if (existingEntry && !existingEntry.id.startsWith("configs/")) {
            return fileId;
        }

        try {
            // Fetch file content
            const body = await this.fetchData<string>(`${this.baseUrl}/${filepath}`, "text");
            const digest = generateDigest(body);

            // Process content based on file type
            const processorConfig = {
                ...config,
                data: {
                    repo: this.repo,
                    ref: this.ref,
                    parent_directory: filepath.includes("/") ? filepath.split("/").slice(0, -1).join("/") : "",
                },
            };

            let rendered: RenderedContent;
            // Use custom processors (may include GitHub-specific transformations for markdown)
            // The processor proxy will handle missing processors appropriately
            rendered = await createProcessors(processors)[extension](body, processorConfig, filepath);

            // Get commit date if requested
            let lastCommit: string | null = null;
            if (getDates) {
                try {
                    const commits = await this.fetchData<GitHubCommit[]>(
                        `https://api.github.com/repos/${this.org}/${this.repo}/commits?path=${filepath}&ref=${this.ref}`,
                        "json",
                    );
                    lastCommit = commits[0]?.commit.committer.date || null;
                } catch {
                    // Ignore commit date errors
                }
            }

            // Store processed file
            store.set({
                id: fileId,
                data: {
                    id: fileId,
                    extension,
                    org: this.org,
                    repo: this.repo,
                    ref: this.ref,
                    lastCommit,
                    ...additionalData,
                },
                body,
                rendered,
                digest,
            });

            return fileId;
        } catch (error) {
            logger.error(`Error processing file ${filepath}: ${error.message}`);
            return undefined;
        }
    }

    /**
     * Process multiple files from a GitHub tree
     */
    async processFilesFromTree({
        path,
        processors,
        getDates = false,
        additionalData = {},
        concurrency = 5,
    }: {
        path: string | RegExp;
        processors: Processors;
        getDates?: boolean;
        additionalData?: Record<string, unknown>;
        concurrency?: number;
    }): Promise<void> {
        const { logger, meta, store } = this.context;
        logger.debug(`Loading files from ${this.org}/${this.repo}@${this.ref}`);

        // Get tree data from GitHub
        const treeData = await this.fetchData<{ sha: string; tree: GithubTreeLeaf[] }>(
            `https://api.github.com/repos/${this.org}/${this.repo}/git/trees/${this.ref}?recursive=1`,
            "json",
        );

        // Skip if tree hasn't changed (scope by repo and ref)
        const treeShaKey = `tree:${this.repo}:${this.ref}`;
        const lastTreeSha = meta.get(treeShaKey);
        if (lastTreeSha === treeData.sha) {
            return;
        }

        // Process matching files
        const matchingLeaves = treeData.tree.filter((leaf) => {
            if (leaf.type === "tree" || leaf.path?.includes(".github/")) return false;
            return typeof path === "string" ? leaf.path?.includes(path) : path.test(leaf.path ?? "");
        });

        const processedIds = new Set<string>();
        const ids = await parallelLimit(
            matchingLeaves.map((leaf) => () => {
                logger.debug(`Processing ${leaf.path}`);
                return this.processFile(leaf.path, processors, getDates, additionalData);
            }),
            concurrency,
        );
        for (const id of ids) {
            if (id) processedIds.add(id);
        }

        // Clean up deleted files (only for this repo/ref)
        const repoRefPrefix = `${this.repo}/${this.ref}/`;
        for (const id of store.keys()) {
            if (id.startsWith(repoRefPrefix) && !processedIds.has(id)) {
                logger.debug(`Removing deleted file: ${id}`);
                store.delete(id);
            }
        }

        // Update tree SHA
        meta.set(treeShaKey, treeData.sha);
        logger.debug("GitHub file processing completed successfully");
    }
}

// ========================================
// AWS S3 HELPERS
// ========================================

/**
 * Helper function to get keys with prefixes from AWS S3
 */
async function getS3KeysWithPrefixes(client: S3Client, prefixes: string[]) {
    const keys: S3Object[] = [];
    const commonPrefixes: S3Prefix[] = [];

    for (const prefix of prefixes) {
        let continuationToken: string | undefined = undefined;

        do {
            const command = new ListObjectsV2Command({
                Bucket: "nf-core-awsmegatests",
                Prefix: prefix,
                ContinuationToken: continuationToken,
                Delimiter: "/",
            });

            try {
                const response = await client.send(command);

                if (response.Contents) {
                    keys.push(...(response.Contents as unknown as S3Object[]));
                }

                if (response.CommonPrefixes) {
                    commonPrefixes.push(...(response.CommonPrefixes as unknown as S3Prefix[]));
                }

                continuationToken = (response as any).NextContinuationToken;
            } catch (error) {
                console.error(`Error retrieving keys for prefix ${prefix}:`, error);
                return { keys: [], commonPrefixes: [] };
            }
        } while (continuationToken);
    }

    return { keys, commonPrefixes };
}

// ========================================
// LOADER FUNCTIONS
// ========================================

/**
 * Parses nf-core config files and extracts profile metadata
 */
function createConfigProcessor(): Processors {
    return {
        config: async (text: string): Promise<RenderedContent> => {
            const NFConfig = {
                config_profile_description: "",
                config_profile_contact: "",
                config_profile_url: "",
                executor: "",
            };

            const descMatch = text.match(/config_profile_description\s*=\s*("+|')([\s\S]*?)\1/);
            NFConfig.config_profile_description = descMatch?.[2] ?? "";

            const contactMatch = text.match(/config_profile_contact\s*=\s*(.*)/);
            NFConfig.config_profile_contact = contactMatch?.[1] ?? "";

            const urlMatch = text.match(/config_profile_url\s*=\s*(.*)/);
            NFConfig.config_profile_url = urlMatch?.[1] ?? "";

            const execMatch = text.match(/executor\s*=\s*(.*)/);
            NFConfig.executor = execMatch?.[1] ?? "";
            if (NFConfig.executor === "") {
                const execNameMatch = text.match(/executor\s*{\s*name\s*=\s*(.*)\s*/);
                NFConfig.executor = execNameMatch?.[1] ?? "";
            }
            if (NFConfig.executor?.includes("?") && NFConfig.executor?.includes(":")) {
                const match = NFConfig.executor.match(/{.*\?(.*):(.*)}/);
                NFConfig.executor = match
                    ? match
                          .slice(1, 3)
                          .map((item) => item.trim().replace(/'/g, ""))
                          .join(", ")
                    : "";
                NFConfig.executor = NFConfig.executor.replace(/\/\/.*$/, "");
            }
            Object.keys(NFConfig).forEach((key) => {
                NFConfig[key] = NFConfig[key]?.replace(/'|"/g, "");
            });
            NFConfig.config_profile_description = NFConfig.config_profile_description.replace(
                " provided by nf-core/configs",
                "",
            );

            return { html: text, metadata: NFConfig };
        },
    };
}

/**
 * Config Loader - Loads nf-core institutional config files
 */
export function configLoader(): Loader {
    return {
        name: "config-loader",
        load: async (context: LoaderContext) => {
            console.log("Starting config loader");
            getCurrentRateLimitRemaining();

            const fetcher = new GitHubContentFetcher("nf-core", "configs", "master", context);

            const processors = {
                md: async (text: string, config: AstroConfig): Promise<RenderedContent> => {
                    const processor = await createMarkdownProcessor(config.markdown);
                    try {
                        const { code: html, metadata } = await processor.render(text);
                        return { html, metadata };
                    } catch (error) {
                        console.error(`Error processing markdown: ${error.message}`);
                        throw error;
                    }
                },
                ...createConfigProcessor(),
            };

            await fetcher.processFilesFromTree({
                path: /\.(md|mdx|config)$/,
                processors,
                getDates: true,
                concurrency: 10,
            });

            console.log("Config loader completed successfully");
        },
    };
}

/**
 * GitHub File Loader - Loads files from GitHub repositories
 */
export function githubFileLoader({ org, repo, ref, processors, path, getDates = false }: GithubLoaderConfig): Loader {
    return {
        name: "github-file-loader",
        load: async (context: LoaderContext) => {
            const fetcher = new GitHubContentFetcher(org, repo, ref, context);
            await fetcher.processFilesFromTree({ path, processors, getDates });
        },
    };
}

/**
 * Pipeline Loader - Loads nf-core pipeline documentation and schemas
 */
export function pipelineLoader(pipelines_json: PipelineJson): Loader {
    return {
        name: "pipeline-loader",
        load: async (context: LoaderContext) => {
            console.log("Starting pipeline loader");
            getCurrentRateLimitRemaining();

            // Flatten all pipeline+release combinations and process concurrently
            const allReleases = pipelines_json.remote_workflows.flatMap((pipeline) =>
                pipeline.releases.map((release) => ({ pipeline, release })),
            );

            await parallelLimit(
                allReleases.map(({ pipeline, release }) => async () => {
                    const fetcher = new GitHubContentFetcher("nf-core", pipeline.name, release.tag_name, context);

                    const filesToProcess = [...release.doc_files, "README.md"];
                    if (release.has_schema) {
                        filesToProcess.push("nextflow_schema.json");
                    }

                    const filePattern = new RegExp(`^(${filesToProcess.map((f) => f.replace(".", "\\.")).join("|")})$`);

                    const processors = {
                        ...createGitHubMarkdownProcessor(context.renderMarkdown, {
                            repo: pipeline.name,
                            ref: release.tag_name,
                        }),
                        ...createSchemaProcessor(context.renderMarkdown),
                    };

                    const metadata = {
                        name: pipeline.name,
                        archived: pipeline.archived,
                        releases: pipeline.releases,
                        description: pipeline.description,
                        topics: pipeline.topics,
                        stargazers_count: pipeline.stargazers_count,
                    };

                    try {
                        await fetcher.processFilesFromTree({
                            path: filePattern,
                            processors,
                            getDates: false,
                            additionalData: metadata,
                        });
                    } catch (error) {
                        context.logger.error(`Error processing ${pipeline.name}@${release.tag_name}: ${error.message}`);
                    }
                }),
                10,
            );

            console.log("Pipeline loader completed successfully");
        },
    };
}

/**
 * Release Notes Loader - Loads GitHub release notes for pipelines
 */
export function releaseLoader(pipelines_json: PipelineJson): Loader {
    return {
        name: "release-notes-loader",
        load: async (context: LoaderContext) => {
            console.log("Starting release notes loader");
            getCurrentRateLimitRemaining();

            const { store, generateDigest, logger, renderMarkdown } = context;

            for (const pipeline of pipelines_json.remote_workflows) {
                logger.debug(`Loading release notes for ${pipeline.name}`);

                try {
                    // Fetch all releases for this pipeline
                    const releases = await octokit.request(`GET /repos/nf-core/${pipeline.name}/releases`);

                    for (const release of releases.data) {
                        const releaseId = `${pipeline.name}/${release.tag_name}/release-notes`;

                        // Skip if unchanged
                        const existingEntry = store.get(releaseId);
                        if (existingEntry && existingEntry.data.tag_name === release.tag_name) {
                            logger.debug(`Release notes for ${pipeline.name}@${release.tag_name} unchanged, skipping`);
                            continue;
                        }

                        const body = release.body || "";
                        const digest = generateDigest(body);

                        try {
                            const { html, metadata } = await renderMarkdown(body);

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
                                    author: release.author
                                        ? {
                                              login: release.author.login,
                                              avatar_url: release.author.avatar_url,
                                              html_url: release.author.html_url,
                                          }
                                        : null,
                                    archived: pipeline.archived,
                                    description: pipeline.description,
                                    topics: pipeline.topics,
                                    stargazers_count: pipeline.stargazers_count,
                                },
                                body,
                                rendered: { html, metadata },
                                digest,
                            });

                            logger.debug(
                                `Successfully processed release notes for ${pipeline.name}@${release.tag_name}`,
                            );
                        } catch (error) {
                            logger.error(
                                `Error processing release notes for ${pipeline.name}@${release.tag_name}: ${error.message}`,
                            );
                        }
                    }
                } catch (error) {
                    logger.error(`Error fetching releases for ${pipeline.name}: ${error.message}`);
                }
            }

            logger.info("GitHub release notes processing completed successfully");
        },
    };
}

/**
 * AWS Results Loader - Fetches pipeline results from AWS S3
 */
export function awsResultsLoader(pipelines_json: PipelineJson): Loader {
    return {
        name: "aws-results-loader",
        load: async (context: LoaderContext) => {
            const { logger, store } = context;
            logger.info("Starting AWS S3 results loader");

            try {
                // Create S3 client (anonymous access)
                const client = new S3Client({
                    region: "eu-west-1",
                    signer: { sign: async (request) => request },
                    credentials: {
                        accessKeyId: "",
                        secretAccessKey: "",
                    },
                });

                const progressBar = new ProgressBar("Processing pipeline: :pipeline [:bar] :percent :etas", {
                    total: pipelines_json.remote_workflows.length,
                });

                // Process each pipeline and its releases
                for (const pipeline of pipelines_json.remote_workflows) {
                    progressBar.tick({ pipeline: pipeline.name });

                    for (const release of pipeline.releases.filter((rel) => rel.tag_name !== "dev")) {
                        const results_path = `results-${release.tag_sha}`;
                        const version = release.tag_name;
                        const pipelineName = pipeline.name;
                        const entryId = `${pipelineName}/${version}/${results_path}`;
                        const prefix = `${pipelineName}/${results_path}/`;

                        // Fetch data from S3
                        const { keys: bucketContents, commonPrefixes } = await getS3KeysWithPrefixes(client, [prefix]);

                        // Skip if no content found
                        if (bucketContents.length === 0 && commonPrefixes.length === 0) {
                            logger.info(`No data found for ${prefix}`);
                            continue;
                        }

                        // Store the entry
                        store.set({
                            id: entryId,
                            data: {
                                pipeline: pipelineName,
                                version: version,
                                results_path: results_path,
                                content: bucketContents,
                                commonPrefixes: commonPrefixes,
                            },
                        });
                    }
                }

                logger.info("AWS results loader completed successfully");
            } catch (error) {
                logger.error(`Error in AWS results loader: ${error.message}`);
            }
        },
    };
}
