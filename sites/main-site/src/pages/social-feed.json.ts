import { getCollection } from "astro:content";
import type { APIContext } from "astro";

// Normalised feed of recently-published / upcoming content, consumed by
// bin/social-post.js to auto-post to Bluesky and Mastodon. Content collections
// are the single source of truth here so we never re-parse frontmatter or
// re-derive URLs elsewhere.

const RECENT_WINDOW_DAYS = 30; // how far back blog posts / advisories stay in the feed
const EVENT_PAST_GRACE_DAYS = 2; // keep just-passed events briefly so same-day posts still fire

type FeedType = "blog" | "advisory" | "event";

interface FeedItem {
    type: FeedType;
    id: string;
    url: string;
    title: string;
    subtitle: string;
    image: string;
    startDate?: string; // events only — the date the auto-poster keys "tomorrow/today" off
}

export async function GET(context: APIContext) {
    const site = context.site!;
    const now = Date.now();
    const recentCutoff = now - RECENT_WINDOW_DAYS * 864e5;
    const eventCutoff = now - EVENT_PAST_GRACE_DAYS * 864e5;

    const slug = (id: string) => id.replace(/\.[^/.]+$/, "");
    const abs = (path: string) => new URL(path, site).href;
    const image = (img?: string) => (img ? abs(img) : abs("/social_img.png"));

    // headerImage is absent on advisories and optional on events; the image()
    // fallback covers both, so one builder serves all three content types.
    const make = (type: FeedType, entry: { id: string; data: any }, extra: Partial<FeedItem> = {}): FeedItem => ({
        type,
        id: `${type}:${slug(entry.id)}`,
        url: abs(`/${type === "advisory" ? "advisories" : type}/${slug(entry.id)}`),
        title: entry.data.title,
        subtitle: entry.data.subtitle,
        image: image(entry.data.headerImage),
        ...extra,
    });

    const blog = await getCollection("blog", ({ data }) => !data.draft && data.pubDate.getTime() <= now);
    const advisories = await getCollection("advisories", ({ data }) => data.publishedDate.getTime() <= now);
    const events = await getCollection("events", ({ data }) => data.start.getTime() >= eventCutoff);

    const items: FeedItem[] = [
        ...blog.filter((p) => p.data.pubDate.getTime() >= recentCutoff).map((p) => make("blog", p)),
        ...advisories.filter((a) => a.data.publishedDate.getTime() >= recentCutoff).map((a) => make("advisory", a)),
        ...events.map((e) => make("event", e, { startDate: e.data.startDate })),
    ];

    return new Response(JSON.stringify(items), {
        headers: { "Content-Type": "application/json" },
    });
}
