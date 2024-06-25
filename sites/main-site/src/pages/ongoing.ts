import type { APIRoute } from 'astro';
import { getCollection, type CollectionEntry } from 'astro:content';
import { EventIsOngoing, newBlogPost } from '@components/store';

export const prerender = false;
export const GET: APIRoute = async ({}) => {
    const now = new Date(Date.now());
    let newPost: CollectionEntry<'blog'> | undefined;
    newPost = await getCollection('blog').then((posts) => {
        return posts.find(
            (post) =>
                now > new Date(post.data.pubDate) &&
                new Date(post.data.pubDate) > new Date(now.getTime() - 1 * 24 * 60 * 60 * 1000),
        );
    });
    if (newPost) {
        newBlogPost.set(true);
    } else {
        newBlogPost.set(false);
    }

    let events: CollectionEntry<'events'>[] = await getCollection('events');
    events = events
        .map((event) => {
            event.data.start = new Date(event.data.startDate + 'T' + event.data.startTime);
            event.data.end = new Date(event.data.endDate + 'T' + event.data.endTime);
            return event;
        })
        .sort((a, b) => {
            return a.data.start - b.data.start;
        });
    events = events.filter((event) => {
        return event.slug.split('/').length === 2;
    });
    const ongoingEvents = (events = events.filter((event) => {
        return event.data.start.getTime() < now.getTime() && event.data.end.getTime() > now.getTime();
    }));
    if (ongoingEvents.length > 0) {
        EventIsOngoing.set(true);
    } else {
        EventIsOngoing.set(false);
    }

    return new Response(JSON.stringify({ newPost: newPost, ongoingEvents: ongoingEvents }), {
        status: 200,
        headers: {
            'Content-Type': 'application/json',
        },
    });
};
