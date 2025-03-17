import rss from '@astrojs/rss';
import { getCollection } from 'astro:content';

export async function GET(context) {
  let blogPosts = await getCollection('blog', ({ data }) => {
    return data.pubDate < new Date();
  });
  blogPosts = blogPosts.sort((a, b) => new Date(b.data.pubDate).getTime() - new Date(a.data.pubDate).getTime());

  return rss({
    title: 'nf-core blog',
    description: 'News and updates from the nf-core community.',
    site: context.site,
    stylesheet: '/rss/styles.xsl',
    items: blogPosts.map((post) => ({
      title: post.data.title,
      description: post.data.subtitle,
      link: '/blog/' + post.id.replace(/\.[^/.]+$/, ''),
      pubDate: post.data.pubDate,
      content: post.body,
    })),
  });
}
