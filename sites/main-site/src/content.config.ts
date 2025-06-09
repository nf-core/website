import { z, defineCollection } from 'astro:content';
import { glob } from 'astro/loaders';
import type { AstroConfig } from 'astro';
import { githubFileLoader } from '@utils/loaders';
import { octokit } from "@components/octokit.js";

const teams = await octokit.request('GET /orgs/{org}/teams', {
    org: 'nf-core',
    headers: {
        'X-GitHub-Api-Version': '2022-11-28'
    }
});

const events = defineCollection({
    loader: glob({ pattern: '**/[^_]*.{md,mdx}', base: './src/content/events' }),
    schema: z
        .object({
            title: z.string(),
            subtitle: z.string(),
            shortTitle: z.string().optional(),
            type: z.enum(['bytesize', 'talk', 'hackathon', 'training']),
            startDate: z.string().refine((s) => /^(\d{4}-\d{2}-\d{2})$/.test(s), {
                message: 'startDate must be in the format YYYY-MM-DD',
            }),
            // check that it contains a time offset
            startTime: z.string().refine((s) => /^(\d{2}:\d{2})([+-]\d{2}:\d{2})$/.test(s), {
                message: 'startTime must be in the format HH:MM+|-HH:MM where the +/-HH:MM is the UTC offset',
            }),
            endDate: z.string().refine((s) => /^(\d{4}-\d{2}-\d{2})$/.test(s), {
                message: 'endDate must be in the format YYYY-MM-DD',
            }),
            endTime: z.string().refine((s) => /^(\d{2}:\d{2})([+-]\d{2}:\d{2})$/.test(s), {
                message: 'endTime must be in the format HH:MM+|-HH:MM where the +/-HH:MM is the UTC offset',
            }),
            announcement: z
                .object({
                    text: z.string().optional(),
                    start: z.date().optional(),
                    end: z.date().optional(),
                })
                .optional(),
            locations: z
                .array(
                    z.object({
                        name: z.string().optional(),
                        links: z.string().url().or(z.string().startsWith('#')).or(z.array(z.string().url())).optional(),
                        geoCoordinates: z.array(z.number(), z.number()).optional(),
                        address: z.string().optional(),
                        country: z.string().optional(),
                        city: z.string().optional(),
                    }),
                )
                .optional(),
            links: z.array(z.string().url()).optional(),
            start: z.date().optional(),
            end: z.date().optional(),
            duration: z.string().optional(),
            embedAt: z.string().optional(),
            importTypeform: z.boolean().optional(),
            hackathonProjectListModals: z.string().optional(),
            youtubeEmbed: z.array(z.string().url()).optional().or(z.string().url()).optional(),
            hideExportButton: z.boolean().optional(),
        })
        .refine((data) => {
            // create start and end date objects
            try {
                data.start = data.start ?? new Date(data.startDate + 'T' + data.startTime);
                data.end = data.end ?? new Date(data.endDate + 'T' + data.endTime);
            } catch (e) {
                throw new Error('startDate and startTime must be in the format YYYY-MM-DD and HH:MM+|-HH:MM');
            }
            // check that start date is before end date
            if (data.start.getTime() > data.end.getTime()) {
                throw new Error(`startDate ${data.start} must be before endDate ${data.end}`);
            }

            // check that announcement.start is before announcement.end
            if (data.announcement?.start && data.announcement.end) {
                if (data.announcement.start.getTime() > data.announcement.end.getTime()) {
                    throw new Error('announcement.start must be before announcement.end');
                }
            }
            // check that announcement.start is set if announcement.text is
            if (data.announcement?.text && !data.announcement.start && !data.announcement.end) {
                throw new Error('announcement.start and announcement.end must be set if announcement.text is');
            }
            // check that locations country is set if locations city is set
            if (data.locations?.[0]?.city && !data.locations?.[0]?.country) {
                throw new Error('locations.country must be set if locations.city is');
            }
            // Return true if the validation should pass
            return true;
        })
        .transform((data) => {
            return {
                ...data,
                start: data.start!,  // Assert dates are set after refinement
                end: data.end!
            };
        }),
});

const about = defineCollection({
    loader: glob({ pattern: '**/[^_]*.{md,mdx}', base: './src/content/about' }),
    schema: z.object({
        title: z.string(),
        description: z.string(),
        md_github_url: z.string().url().optional(),
        minHeadingDepth: z.number().optional(),
        maxHeadingDepth: z.number().optional(),
    }),
});

const blog = defineCollection({
    loader: glob({ pattern: '**/[^_]*.{md,mdx}', base: './src/content/blog' }),
    schema: z
        .object({
            title: z.string(),
            subtitle: z.string(),
            shortTitle: z.string().optional(),
            headerImage: z.string().url().optional().or(z.string().startsWith('/assets/images/blog/')).optional(),
            headerImageAlt: z.string().optional(),
            headerImageDim: z.array(z.number(), z.number()).optional(),
            label: z.array(z.string()),
            pubDate: z.date(),
            authors: z.array(z.string()),
            draft: z.boolean().optional(),
            embedHeaderImage: z.boolean().optional(),
            announcement: z
                .object({
                    text: z.string().optional(),
                    start: z.date().optional(),
                    end: z.date().optional(),
                })
                .optional(),
            maxHeadingDepth: z.number().optional(),
        })
        .refine((data) => {
            // Check if headerImage is present but headerImageAlt is not
            if (data.headerImage && !data.headerImageAlt) {
                throw new Error('Please provide alt text for your `headerImage` in `headerImageAlt`.');
            }
            // Check if headerImageDim is present but headerImage is not present or does not start with /assets/
            if (data.headerImageDim && (!data.headerImage || !data.headerImage.startsWith('/assets/'))) {
                throw new Error(
                    'Please provide a `headerImage` that starts with `/assets/` if you are providing `headerImageDim`.',
                );
            }
            // check that announcement.start is before announcement.end
            if (data.announcement?.start && data.announcement.end) {
                if (data.announcement.start.getTime() > data.announcement.end.getTime()) {
                    throw new Error('`announcement.start` must be before `announcement.end`');
                }
            }
            // check that announcement.start is set if announcement.text is
            if (data.announcement?.text && !data.announcement.start && !data.announcement.end) {
                throw new Error('`announcement.start` and `announcement.end` must be set if `announcement.text` is');
            }
            // Return true if the validation should pass
            return true;
        }),
});

const specialInterestGroups = defineCollection({
    loader: glob({ pattern: '**/[^_]*.{md,mdx}', base: './src/content/special-interest-groups' }),
    schema: z
        .object({
            title: z.string(),
            subtitle: z.string(),
            groupName: z.string(),
            // for index.md pages also require lead and pipelines
            leads: z
                .array(z.string())
                .or(z.array(z.record(z.string())))
                .optional(),
            pipelines: z
                .array(z.string())
                .optional()
                .transform((data) => {
                    // sort the pipelines by name
                    return data?.sort();
                }),
        })
        .refine((data) => {
            if (data?.leads && !data.pipelines) {
                throw new Error('`pipelines` must be set if `leads` is');
            }
            return true;
        }),
});

const hackathonProjects = defineCollection({
    loader: glob({ pattern: '**/[^_]*.{md,mdx}', base: './src/content/hackathon-projects' }),
    schema: z
        .object({
            title: z.string(),
            category: z.enum(['pipelines', 'components', 'tooling', 'community']),
            leaders: z.record(z.object({
                name: z.string(),
                slack: z.string().url().refine((url) => {
                    if (url && !url.startsWith('https://nfcore.slack.com/')) {
                        throw new Error('leaders `slack` must be a nf-core or Nextflow slack URL: '+url);
                    }
                    return true;
                }).optional()
            })),
            color: z
                .string()
                .refine((data) => {
                    if (data && !data.startsWith('#') && !data.startsWith("'#")) {
                        throw new Error('`color` must start with "#"');
                    }
                    return true;
                })
                .optional(),
            intro_video: z.string().optional(),
            image: z.string().optional(),
            image_alt: z.string().optional(),
            slack: z.string().url().refine((url) => {
                    if (!url.startsWith('https://nfcore.slack.com/')) {
                        throw new Error('`slack` must be a nf-core slack URL: '+url);
                    }
                    return true;
                }),
        })
        .refine((data) => {
            // Check if headerImage is present but headerImageAlt is not
            if (data.image && !data.image_alt) {
                throw new Error('Please provide alt text for your `image` in `image_alt`.');
            }
            // Return true if the validation should pass
            return true;
        }),
});

const stickers = defineCollection({
    loader: githubFileLoader({
        org: 'nf-core',
        repo: 'logos',
        ref: 'master',
        path: /^hexagon-stickers\/.+\.png$/,
        getDates: false,
        processors: {
            png: async (text: string, config: AstroConfig, filepath?: string) => {
                if (!filepath) {
                    console.log("filepath is undefined");
                    return {
                        html: '',
                        metadata: {}
                    };
                }
                // skip template and preview images
                if (filepath.includes('-TEMPLATE-') || filepath.includes('-preview')) {
                    console.log("skipping template or preview image", filepath);
                    return {
                        html: '',
                        metadata: {}
                    };
                }
                // Extract category and name from the file path
                const pathMatch = filepath.match(/^hexagon-stickers\/([^/]+)(?:\/([^/]+))?\/[^/]+\.png$/);
                if (!pathMatch) {
                    return {
                        html: '',
                        metadata: {
                            name: 'unknown',
                            category: 'other',
                            imageUrl: '',
                            alt: 'Unknown sticker',
                            link: 'https://nf-co.re',
                            description: 'Unknown sticker'
                        }
                    };
                }

                const [, category, subdir] = pathMatch;
                const name = subdir || category;
                const imageUrl = `https://raw.githubusercontent.com/nf-core/logos/master/${filepath}`;
                let description = subdir ? `${name} pipeline sticker` : '';
                if (category === 'pipelines') {
                    const response = await octokit.request('GET /repos/{owner}/{repo}', {
                        owner: 'nf-core',
                        repo: name,
                        headers: {
                            'X-GitHub-Api-Version': '2022-11-28'
                        }
                    });
                    description = response.data.description || description;
                }
                if (category === 'teams') {
                    description = teams.data.find((team) => team.name === filepath.split('/').pop()?.replace("nf-core-", "").replace(".png", ""))?.description || description;
                }

                return {
                    html: '',
                    metadata: {
                        name,
                        category,
                        imageUrl,
                        alt: `${name} sticker`,
                        link: subdir ? `https://github.com/nf-core/${name}` : 'https://nf-co.re',
                        description
                    }
                };
            }
        }
    }),
    schema: z.object({
        name: z.string(),
        category: z.string(),
        imageUrl: z.string(),
        alt: z.string(),
        link: z.string().url(),
        description: z.string()
    })
});

export const collections = {
    events,
    about,
    blog,
    'special-interest-groups': specialInterestGroups,
    'hackathon-projects': hackathonProjects,
    stickers
};
