import rss from '@astrojs/rss';
import { getCollection } from 'astro:content';

export async function GET(context) {
  const advisories = await getCollection('advisories');
  return rss({
    title: 'nf-core: Advisories',
    description: 'nf-core advisories and security announcements.',
    site: context.site,
    stylesheet: '/rss/styles.xsl',
    items: advisories.map((advisory) => ({
      title: advisory.data.title,
      description: advisory.data.subtitle,
      link: '/advisories/' + advisory.id.replace(/\.[^/.]+$/, ''),
      pubDate: advisory.data.published_date,
    })),
  });
}
