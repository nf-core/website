import admonitionsPlugin from '../../bin/remark-admonitions';
import { mermaid } from '../../bin/remark-mermaid.ts';
import { rehypeCheckboxParser } from '../../bin/rehype-checkbox-parser.ts';
import { rehypeHeadingNumbers } from '../../bin/rehype-heading-numbers.ts';
import mdx from '@astrojs/mdx';
import netlify from '@astrojs/netlify';
import partytown from '@astrojs/partytown';
import sitemap from '@astrojs/sitemap';
import svelte from '@astrojs/svelte';
import yaml from '@rollup/plugin-yaml';
import { defineConfig, envField } from 'astro/config';
import { FontaineTransform } from 'fontaine';
import { h } from 'hastscript';
import addClasses from 'rehype-class-names';
import rehypeAutolinkHeadings from 'rehype-autolink-headings';
import rehypeKatex from 'rehype-katex';
import rehypePrettyCode from 'rehype-pretty-code';
import rehypeSlug from 'rehype-slug';
import urls from 'rehype-urls';
import rehypeWrap from 'rehype-wrap-all';
import remarkDirective from 'remark-directive';
import emoji from 'remark-emoji';
import remarkGfm from 'remark-gfm';
import remarkMath from 'remark-math';
import pipelines_json from './public/pipelines.json';
import markdownIntegration from '@astropub/md';
import icon from 'astro-icon';

let pipelineRedirects = {};
pipelines_json.remote_workflows.map((pipeline) => {
    pipelineRedirects[`/${pipeline.name}/:version/*`] =
        `https://nf-core-pipelines.netlify.app/${pipeline.name}/:version/:splat 200!`;
    pipelineRedirects[`/${pipeline.name}/`] = `https://nf-core-pipelines.netlify.app/${pipeline.name} 200!`;
});
// https://astro.build/config
export default defineConfig({
	site: 'https://nf-co.re/',
	output: 'static',
	adapter: netlify(),
	prefetch: false,
	redirects: {
		...pipelineRedirects
	},
	env: {
		schema: {
			GITHUB_TOKEN: envField.string({
				context: 'server',  // Keep as server-side only for security
                access: 'secret',
                optional: false
			})
		}
	},
	integrations: [
		svelte(),
		icon({
			include: {
				// only include a subset of icons
				'file-icons': ['nextflow'],
				logos: [
					'twitter',
					'mastodon-icon',
					'slack-icon',
					'aws',
					'microsoft-azure',
					'github-actions',
					'youtube-icon',
					'linkedin'
				],
				fa: ['github'],
				'fa-brands': ['github'],
				mdi: [
					'aws',
					'slack',
					'youtube',
					'cloud-outline',
					'timeline-check-outline',
					'book-information-variant',
					'package-variant',
					'progress-check'
				],
				octicon: [
					'chevron-right-16',
					'git-pull-request-16',
					'law-16',
					'link-external-16',
					'mortar-board-16',
					'play-16',
					'table-16',
					'tasklist-16',
					'terminal-16',
					'tools-16'
				],
				'simple-icons': ['bluesky'],
				ri: ['open-source-line'],
				'pepicons-print': ['t-shirt'],
				fluent: ['paint-brush-sparkle-24-filled'],
				bi: ['bag-heart-fill']
			}
		}),
		sitemap(),
		partytown({
			// Adds dataLayer.push as a forwarding-event.
			config: {
				forward: ['dataLayer.push']
			}
		}),
		mdx(),
		markdownIntegration()
	],
	build: {
		inlineStylesheets: 'auto',
		format: 'file',
		assetsPrefix:
			process.env.CONTEXT === 'production'
				? 'https://nf-core-main-site.netlify.app/'
				: process.env.DEPLOY_PRIME_URL
	},
	vite: {
		plugins: [
			yaml(),
			FontaineTransform.vite({
				// avoid flash of unstyled text by interjecting fallback system fonts https://developer.chrome.com/blog/framework-tools-font-fallback/#using-fontaine-library
				fallbacks: [
					'BlinkMacSystemFont',
					'Segoe UI',
					'Helvetica Neue',
					'Arial',
					'Noto Sans'
				],
				resolvePath: (id) => new URL(`./public${id}`, import.meta.url),
				skipFontFaceGeneration: (fallbackName) =>
					fallbackName === 'Font Awesome 6 Pro fallback'
			})
		],
		resolve: {
			preserveSymlinks: true
		},
		envPrefix: ['PUBLIC_', 'GITHUB_'],  // This allows GITHUB_ prefixed env vars
	},
	image: {
		domains: [
			'raw.githubusercontent.com',
			'unsplash.com',
			'netlify.app',
			'nf-co.re',
			'nf-core-main-site.netlify.app'
		],
		service: {
			entrypoint: 'astro/assets/services/sharp'
		}
	},
	markdown: {
		syntaxHighlight: false,
		remarkPlugins: [
			emoji,
			remarkGfm,
			remarkDirective,
			admonitionsPlugin,
			mermaid,
			remarkMath
			// [
			//     remarkDescription,
			//     {
			//         name: 'excerpt',
			//         node: (node, i, parent) => {
			//             // check if parent has a child that is an html comment with the text 'end of excerpt'
			//             if (
			//                 parent?.children?.some(
			//                     (child) =>
			//                         (child.type === 'html' && child.value === '<!-- end of excerpt -->') ||
			//                         (child.type === 'mdxFlowExpression' && child?.value === '/* end of excerpt */'),
			//                 )
			//             ) {
			//                 const sibling = parent?.children[i + 1];

			//                 return (
			//                     (sibling?.type === 'html' && sibling?.value === '<!-- end of excerpt -->') ||
			//                     (sibling?.type === 'mdxFlowExpression' && sibling?.value === '/* end of excerpt */')
			//                 );
			//             } else {
			//                 // return the first paragraph otherwise

			//                 // get the index of the first paragraph
			//                 const firstParagraphIndex = parent?.children.findIndex(
			//                     (child) => child.type === 'paragraph',
			//                 );
			//                 // if the node is the first paragraph, return true
			//                 return i === firstParagraphIndex;
			//             }
			//         },
			//     },
			// ],
		],
		// NOTE: Also update the plugins in `src/components/Markdown.svelte`!
		rehypePlugins: [
			rehypeSlug,
			[
				rehypeAutolinkHeadings,
				{
					behavior: 'append',
					content: h('i.ms-1.fas.fa-link.fa-xs.invisible')
				}
			],
			[
				addClasses,
				{
					table: 'table table-hover table-sm small'
				}
			],
			[
				rehypeWrap,
				{
					selector: 'table',
					wrapper: 'div.table-responsive'
				}
			],
			[
				urls,
				(url) => {
					const regex = /^https:\/\/(raw.)*github/;
					if (!regex.test(url.href) && url.href?.endsWith('.md')) {
						url.href = url.href.replace(/\.md$/, '/');
						url.pathname = url.pathname.replace(/\.md$/, '/');
						url.path = url.path.replace(/\.md$/, '/');
					} else if (!regex.test(url.href) && url.href?.endsWith('.mdx')) {
						url.href = url.href.replace(/\.mdx$/, '/');
						url.pathname = url.pathname.replace(/\.mdx$/, '/');
						url.path = url.path.replace(/\.mdx$/, '/');
					}
				}
			],
			rehypeCheckboxParser,
			rehypeHeadingNumbers,
			[
				rehypePrettyCode,
				{
					defaultLang: 'plaintext',
					keepBackground: true,
					theme: {
						dark: 'github-dark',
						light: 'github-light'
					}
				}
			],
			rehypeKatex
		]
	}
});
