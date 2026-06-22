// Builds the full content-data object for a single newsletter month, including
// cropped image thumbnails. Kept separate from newsletter.ts because it imports
// `astro:assets` (a virtual module only resolvable inside the Astro build) — the
// pure data/date helpers in newsletter.ts must stay importable from plain Node
// (e.g. the Playwright unit tests), which can't resolve the `astro:` scheme.
import { getImage } from "astro:assets";
import {
    getMonthName,
    getBlogPostsForMonth,
    getBlogPostsForPreviousMonths,
    getPipelineReleasesForMonth,
    getNewPipelines,
    getEventsForMonth,
    getFutureEvents,
    getAdvisoriesForMonth,
    getAdvisoriesForPreviousMonths,
    getProposalsForMonth,
} from "./newsletter";
import type { PipelineWorkflow, RawProposal } from "./newsletter";

export async function getNewsletterContentData(
    getCollectionFn: (name: string) => Promise<any[]>,
    pipelines: PipelineWorkflow[],
    year: number,
    month: number,
    allProposals: RawProposal[],
    images: Record<string, () => Promise<any>>,
) {
    const monthName = getMonthName(month);
    const baseUrl = "https://nf-co.re";

    // Content month (previous calendar month — newsletter looks back)
    const contentMonth = month === 1 ? 12 : month - 1;
    const contentYear = month === 1 ? year - 1 : year;

    let blogPosts = await getCollectionFn("blog");
    blogPosts = blogPosts.filter((post: any) => new Date(post.data.pubDate) < new Date());

    let events = await getCollectionFn("events");
    events = events.filter((e: any) => e.id.split("/").length === 2);

    const thisMonthBlogPosts = getBlogPostsForMonth(blogPosts, contentYear, contentMonth);
    const olderBlogPosts = getBlogPostsForPreviousMonths(blogPosts, contentYear, contentMonth, 2);
    const releases = getPipelineReleasesForMonth(pipelines, contentYear, contentMonth);
    const firstReleases = releases.filter((r) => r.isFirstRelease);
    const otherReleases = releases.filter((r) => !r.isFirstRelease);
    const newPipelines = getNewPipelines(pipelines, contentYear, contentMonth);
    const recentEvents = getEventsForMonth(events, contentYear, contentMonth);
    const allAdvisories = await getCollectionFn("advisories");
    const monthAdvisories = getAdvisoriesForMonth(allAdvisories, contentYear, contentMonth);
    const olderAdvisories = getAdvisoriesForPreviousMonths(allAdvisories, contentYear, contentMonth, 2);
    const proposals = getProposalsForMonth(allProposals, contentYear, contentMonth);
    const pipelineProposals = proposals.filter((p) => p.category === "pipeline");
    const otherProposals = proposals.filter((p) => p.category === "other");

    const upcomingEvents = [...getEventsForMonth(events, year, month), ...getFutureEvents(events, year, month, 1)];

    // Build email preview summary (~90 chars for Gmail)
    const previewParts: string[] = [];
    if (thisMonthBlogPosts.length > 0)
        previewParts.push(`${thisMonthBlogPosts.length} blog post${thisMonthBlogPosts.length > 1 ? "s" : ""}`);
    if (monthAdvisories.length > 0)
        previewParts.push(`${monthAdvisories.length} advisor${monthAdvisories.length > 1 ? "ies" : "y"}`);
    if (releases.length > 0) previewParts.push(`${releases.length} pipeline release${releases.length > 1 ? "s" : ""}`);
    if (upcomingEvents.length > 0)
        previewParts.push(`${upcomingEvents.length} upcoming event${upcomingEvents.length > 1 ? "s" : ""}`);
    const previewText =
        previewParts.length > 0 ? previewParts.join(", ") + "." : `Community news for ${monthName} ${year}.`;
    const newsletterUrl = `${baseUrl}/newsletter/${year}/${String(month).padStart(2, "0")}`;

    // Pre-resolve all image sources to a cropped 16:9 thumbnail. Header images
    // come in any aspect ratio, and email clients ignore `object-fit`, so we crop
    // the actual file here (works for both bundled /assets images and whitelisted
    // remote URLs). Falls back to the un-cropped src if optimisation fails (e.g. a
    // non-whitelisted domain), which still renders — just uncropped.
    async function resolveImageSrc(headerImage: string): Promise<string> {
        let src: any = headerImage;
        if (headerImage.startsWith("/assets/")) {
            const loader = images["/src" + headerImage];
            if (!loader) return headerImage;
            src = ((await loader()) as any).default;
        }
        try {
            // JPEG (not WebP) for broad email-client support, incl. Outlook on Windows.
            // JPG (not WebP) for broad email-client support, incl. Outlook on Windows.
            const thumb = await getImage({
                src,
                width: 480,
                height: 270,
                fit: "cover",
                position: "center",
                format: "jpg",
            });
            return thumb.src;
        } catch {
            return typeof src === "string" ? src : src.src;
        }
    }

    // Resolve sequentially, not via Promise.all: getImage lazily imports and
    // caches the image service on first use, and firing the whole batch
    // concurrently races that init (some calls get an undefined service).
    const blogImageSrcs = new Map<string, any>();
    for (const post of [...thisMonthBlogPosts, ...olderBlogPosts]) {
        blogImageSrcs.set(post.id, await resolveImageSrc((post as any).data.headerImage));
    }
    const eventImageSrcs = new Map<string, any>();
    for (const event of upcomingEvents) {
        if ((event as any).data.headerImage) {
            eventImageSrcs.set(event.id, await resolveImageSrc((event as any).data.headerImage));
        }
    }

    return {
        year,
        month,
        monthName,
        baseUrl,
        newsletterUrl,
        previewText,
        thisMonthBlogPosts,
        olderBlogPosts,
        firstReleases,
        otherReleases,
        newPipelines,
        recentEvents,
        upcomingEvents,
        pipelineProposals,
        otherProposals,
        monthAdvisories,
        olderAdvisories,
        blogImageSrcs,
        eventImageSrcs,
    };
}
