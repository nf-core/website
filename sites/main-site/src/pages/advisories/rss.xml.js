import rss from '@astrojs/rss';
import { getCollection } from 'astro:content';

export async function GET(context) {
  let advisories = await getCollection('advisories');
  advisories = advisories.sort((a, b) => new Date(b.data.publishedDate).getTime() - new Date(a.data.publishedDate).getTime());

  return rss({
    title: 'nf-core: Advisories',
    description: 'nf-core advisories and security announcements.',
    site: context.site,
    stylesheet: '/rss/styles.xsl',
    items: advisories.map((advisory) => ({
      title: advisory.data.title,
      description: advisory.data.subtitle,
      link: '/advisories/' + advisory.id.replace(/\.[^/.]+$/, ''),
      pubDate: advisory.data.publishedDate,
    })),
  });
}
