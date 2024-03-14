import rss from '@astrojs/rss';
import { getCollection } from 'astro:content';

export async function GET(context) {
  const events = await getCollection('events');
  return rss({
    title: 'nf-core: Events',
    description: 'Details of past and future nf-core events.',
    site: context.site,
    stylesheet: '/rss/styles.xsl',
    items: events.map((event) => ({
      title: event.data.title,
      description: event.data.subtitle || '',
      pubDate: new Date(event.data.startDate + 'T' + event.data.startTime),
      link: '/events/' + event.slug,
    })),
  });
}
