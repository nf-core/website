import rss from '@astrojs/rss';
import { getCollection } from 'astro:content';

export async function GET(context) {
  const blogPosts = await getCollection('blog');
  return rss({
    title: 'nf-core: Events',
    description: 'News and updates from the nf-core community.',
    site: context.site,
    stylesheet: '/rss/styles.xsl',
    items: blogPosts.map((post) => ({
      title: post.data.title,
      description: post.data.subtitle,
      link: '/blog/' + post.slug,
      pubDate: post.data.pubDate,
    })),
  });
}
