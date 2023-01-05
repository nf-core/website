import rss from '@astrojs/rss';

const eventImportResult = import.meta.glob('./**/*.md', { eager: true });
const events = Object.values(eventImportResult);

export const get = () =>
  rss({
    title: 'nf-core: Events',
    description: 'Details of past and future nf-core meetups.',
    site: import.meta.env.SITE,
    stylesheet: '/rss/styles.xsl',
    items: events.map((event) => ({
      link: event.url,
      title: event.frontmatter.title,
      // TODO: FIX WITH PROPER FRONTMATTER
      // date('r', $event['start_ts'])
      // pubDate: event.frontmatter.start_ts,
    })),
  });
