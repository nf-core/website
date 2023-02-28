// https://astro.build/config
import prefetch from '@astrojs/prefetch';
import sitemap from '@astrojs/sitemap';
import svelte from '@astrojs/svelte';
import markdownIntegration from '@astropub/md';
import yaml from '@rollup/plugin-yaml';
import { defineConfig } from 'astro/config';
import { h } from 'hastscript';
import addClasses from 'rehype-add-classes';
import rehypeAutolinkHeadings from 'rehype-autolink-headings';
import rehypePrettyCode from 'rehype-pretty-code';
import rehypeSlug from 'rehype-slug';
import urls from 'rehype-urls';
import rehypeWrap from 'rehype-wrap-all';
import emoji from 'remark-emoji';
import remarkGfm from 'remark-gfm';
import { BUNDLED_LANGUAGES } from 'shiki';
import partytown from "@astrojs/partytown";
BUNDLED_LANGUAGES = BUNDLED_LANGUAGES.map(lang => {
  if (lang.id === 'groovy') {
    lang.aliases = ['nextflow', 'nf'];
  }
  return lang;
});

// https://astro.build/config
export default defineConfig({
  site: 'https://nf-co.re/',
  integrations: [svelte(), sitemap(), markdownIntegration(), prefetch(), partytown()],
  vite: {
    plugins: [yaml()],
    ssr: {
      noExternal: ['@popperjs/core', 'bin/cache.js']
    }
  },
  markdown: {
    syntaxHighlight: false,
    remarkPlugins: [emoji, remarkGfm],
    rehypePlugins: [rehypeSlug, [rehypeAutolinkHeadings, {
      behavior: 'append',
      content: h('i.ms-1.fas.fa-link.invisible')
    }], [addClasses, {
      table: 'table table-hover table-sm small'
    }], [rehypeWrap, {
      selector: 'table',
      wrapper: 'div.table-responsive'
    }], [urls, url => {
      if (url.href?.endsWith('.md')) {
        url.href = url.href.replace(/\.md$/, '/');
        url.pathname = url.pathname.replace(/\.md$/, '/');
        url.path = url.path.replace(/\.md$/, '/');
      }
    }], [rehypePrettyCode, {
      langPrefix: 'language-'
    }]
    // [rehypeWrap, { selector: 'pre:has(code.language-bash)', wrapper: 'div.copy-code' }],
    ]
  }
});