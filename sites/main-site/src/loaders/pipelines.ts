import type { Loader } from 'astro/loaders';
import pipelines_json from '@public/pipelines.json';
import { createMarkdownProcessor } from '@astrojs/markdown-remark';

export const pipelinesLoader: Loader = {
    name: 'pipelines',
    load: async ({ store, logger, config, generateDigest }) => {
        const processor = await createMarkdownProcessor(config.markdown);
        // store.clear();
        pipelines_json.remote_workflows.map(async (pipeline) => {
            logger.info(`Loaded pipeline: ${pipeline.name}`);
            if (pipeline.archived) {
                return;
            }
            pipeline.releases.map(async (version) => {
                const paths = ['README.md', ...version.doc_files];
                paths.map(async (path) => {
                    const repo = pipeline.name;
                    const ref = version.tag_name;

                    if (store.get(pipeline.name + '-' + ref + '-' + path)) {
                        logger.info(`Skipping ${pipeline.name}-${ref}-${path}`);
                        return;
                    }
                    logger.info(`fetching file: ${ref}=${path}`);
                    const readmeUrl = `https://raw.githubusercontent.com/nf-core/${repo}/${ref}/${path}`;
                    const response = await fetch(readmeUrl).catch((error) => {
                        logger.error(`Error fetching ${readmeUrl}: ${error.message}`);
                    });
                    let content = await response?.text();
                    if (!content) {
                        return;
                    }
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
                    const renderedContent = await processor.render(content);

                    store.set({
                        id: pipeline.name + '-' + ref + '-' + path,
                        digest: generateDigest(content),
                        body: content,
                        rendered: {
                            html: renderedContent.code,
                            metadata: {
                                headings: renderedContent.metadata.headings,
                            },
                        },
                        data: {
                            name: pipeline.name,
                        },
                    });
                });
            });
        });
    },
};
