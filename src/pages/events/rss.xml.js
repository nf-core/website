import rss, { pagesGlobToRssItems } from '@astrojs/rss';

export const get = () =>
  rss({
    title: 'nf-core: Events',
    description: 'Details of past and future nf-core meetups.',
    site: import.meta.env.SITE,
    stylesheet: '/rss/styles.xsl',
    items: pagesGlobToRssItems(import.meta.glob('./events/*.{md,mdx}')),
  });
