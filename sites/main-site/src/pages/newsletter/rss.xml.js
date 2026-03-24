import rss from "@astrojs/rss";
import { getCollection } from "astro:content";
import pipelines_json from "@public/pipelines.json";
import {
    getNewsletterMonths,
    getBlogPostsForMonth,
    getPipelineReleasesForMonth,
    getEventsForMonth,
    getFutureEvents,
    getMonthName,
    getAdvisoriesForMonth,
} from "@utils/newsletter";

export async function GET(context) {
    let blogPosts = await getCollection("blog");
    blogPosts = blogPosts.filter((post) => new Date(post.data.pubDate) < new Date());

    let events = await getCollection("events");
    events = events.filter((e) => e.id.split("/").length === 2);

    const pipelines = pipelines_json.remote_workflows;
    const advisories = await getCollection("advisories");
    const months = getNewsletterMonths(blogPosts, events, pipelines, advisories);

    const items = months.map(({ year, month }) => {
        const monthName = getMonthName(month);
        const contentMonth = month === 1 ? 12 : month - 1;
        const contentYear = month === 1 ? year - 1 : year;

        // Build a summary description
        const parts = [];
        const blogCount = getBlogPostsForMonth(blogPosts, contentYear, contentMonth).length;
        const releaseCount = getPipelineReleasesForMonth(pipelines, contentYear, contentMonth).length;
        const eventCount = [
            ...getEventsForMonth(events, year, month),
            ...getFutureEvents(events, year, month, 1),
        ].length;
        const advisoryCount = getAdvisoriesForMonth(advisories, contentYear, contentMonth).length;

        if (blogCount > 0) parts.push(`${blogCount} blog post${blogCount > 1 ? "s" : ""}`);
        if (advisoryCount > 0) parts.push(`${advisoryCount} advisor${advisoryCount > 1 ? "ies" : "y"}`);
        if (releaseCount > 0) parts.push(`${releaseCount} pipeline release${releaseCount > 1 ? "s" : ""}`);
        if (eventCount > 0) parts.push(`${eventCount} upcoming event${eventCount > 1 ? "s" : ""}`);

        const description =
            parts.length > 0
                ? `nf-core community news: ${parts.join(", ")}.`
                : `nf-core community newsletter for ${monthName} ${year}.`;

        return {
            title: `nf-core Newsletter - ${monthName} ${year}`,
            description,
            link: `/newsletter/${year}/${String(month).padStart(2, "0")}`,
            pubDate: new Date(year, month - 1, 1),
        };
    });

    return rss({
        title: "nf-core Newsletter",
        description: "Monthly community newsletter from the nf-core project.",
        site: context.site,
        stylesheet: "/rss/styles.xsl",
        items,
    });
}
