import { octokit } from "@components/octokit.js";
import { addMonths } from "date-fns";

// ========================================
// TYPE DEFINITIONS
// ========================================

export interface PipelineWorkflow {
    name: string;
    description: string;
    created_at: string;
    archived: boolean;
    topics: string[];
    stargazers_count: number;
    releases: PipelineRelease[];
}

interface PipelineRelease {
    tag_name: string;
    published_at: string;
    tag_sha: string;
    has_schema: boolean;
}

export interface NewsletterMonth {
    year: number;
    month: number;
}

export interface NewsletterRelease {
    pipelineName: string;
    description: string;
    tagNames: string[];
    publishedAt: string;
    isFirstRelease: boolean;
}

export interface NewPipeline {
    name: string;
    description: string;
    createdAt: string;
    isDevOnly: boolean;
}

export interface RawProposal {
    title: string;
    url: string;
    number: number;
    labels: string[];
    closedAt: string | null;
    stateReason: string | null;
    createdAt: string;
}

export interface Proposal extends RawProposal {
    displayTitle: string;
    category: "pipeline" | "other";
    status: "new" | "accepted";
}

// ========================================
// DATE UTILITIES
// ========================================

export function getMonthRange(year: number, month: number): { start: Date; end: Date } {
    const start = new Date(Date.UTC(year, month - 1, 1));
    const end = new Date(Date.UTC(year, month, 1));
    return { start, end };
}

export function getMonthName(month: number): string {
    return new Date(2000, month - 1, 1).toLocaleString("en-GB", { month: "long" });
}

function isInMonth(date: Date, year: number, month: number): boolean {
    const { start, end } = getMonthRange(year, month);
    return date >= start && date < end;
}

/**
 * The `count` calendar months adjacent to (year, month): `step` of -1 walks
 * backwards (preceding months), +1 walks forwards. Anchored at noon UTC so
 * date-fns month maths stays on the right month regardless of host timezone.
 */
function adjacentMonths(year: number, month: number, count: number, step: 1 | -1): NewsletterMonth[] {
    const anchor = new Date(Date.UTC(year, month - 1, 1, 12));
    return Array.from({ length: count }, (_, i) => {
        const d = addMonths(anchor, step * (i + 1));
        return { year: d.getUTCFullYear(), month: d.getUTCMonth() + 1 };
    });
}

// ========================================
// DATA FILTERING
// ========================================

/**
 * Get all months that have any newsletter content, sorted descending.
 */
export function getNewsletterMonths(
    blogPosts: { data: { pubDate: Date } }[],
    events: { data: { start: Date } }[],
    pipelines: PipelineWorkflow[],
    advisories: { data: { publishedDate: Date } }[] = [],
): NewsletterMonth[] {
    const monthSet = new Set<string>();

    for (const post of blogPosts) {
        const d = new Date(post.data.pubDate);
        monthSet.add(`${d.getUTCFullYear()}-${d.getUTCMonth() + 1}`);
    }

    for (const event of events) {
        const d = new Date(event.data.start);
        monthSet.add(`${d.getUTCFullYear()}-${d.getUTCMonth() + 1}`);
    }

    for (const advisory of advisories) {
        const d = new Date(advisory.data.publishedDate);
        monthSet.add(`${d.getUTCFullYear()}-${d.getUTCMonth() + 1}`);
    }

    for (const pipeline of pipelines) {
        for (const release of pipeline.releases) {
            if (release.tag_name === "dev") continue;
            const d = new Date(release.published_at);
            monthSet.add(`${d.getUTCFullYear()}-${d.getUTCMonth() + 1}`);
        }
        // Also include month when pipeline was created
        const created = new Date(pipeline.created_at);
        monthSet.add(`${created.getUTCFullYear()}-${created.getUTCMonth() + 1}`);
    }

    const months: NewsletterMonth[] = [...monthSet].map((key) => {
        const [year, month] = key.split("-").map(Number);
        return { year, month };
    });

    // Sort descending (newest first)
    months.sort((a, b) => {
        if (a.year !== b.year) return b.year - a.year;
        return b.month - a.month;
    });

    // Filter to past/current months only
    const now = new Date();
    const currentYear = now.getUTCFullYear();
    const currentMonth = now.getUTCMonth() + 1;
    return months.filter((m) => m.year < currentYear || (m.year === currentYear && m.month <= currentMonth));
}

/**
 * Whether the newsletter dated the 1st of (year, month) would render any content.
 *
 * The candidate-month set above is seeded from the months that content lands in,
 * but a newsletter looks *back* at the previous calendar month (blog posts,
 * advisories, releases, new pipelines, recent events, proposals) and *forward* at
 * this + next month (upcoming events). That offset means a seeded month can end up
 * with every section empty — e.g. a lone advisory dated in a month whose preceding
 * month has nothing. This mirrors the section gating in NewsletterLayout so we only
 * keep months where at least one section would actually render, and never publish a
 * page that is just the header and footer.
 */
export function newsletterMonthHasContent(
    blogPosts: { id: string; data: { pubDate: Date } }[],
    events: { id: string; data: { start: Date } }[],
    pipelines: PipelineWorkflow[],
    advisories: { id: string; data: { publishedDate: Date } }[],
    proposals: RawProposal[],
    year: number,
    month: number,
): boolean {
    // Previous calendar month, via the date-fns helper so the year boundary is
    // handled the same way as everywhere else (no manual `month === 1 ? …` maths).
    const [{ year: contentYear, month: contentMonth }] = adjacentMonths(year, month, 1, -1);

    const checks = [
        () => getBlogPostsForMonth(blogPosts, contentYear, contentMonth),
        () => getBlogPostsForPreviousMonths(blogPosts, contentYear, contentMonth, 2),
        () => getAdvisoriesForMonth(advisories, contentYear, contentMonth),
        () => getAdvisoriesForPreviousMonths(advisories, contentYear, contentMonth, 2),
        () => getPipelineReleasesForMonth(pipelines, contentYear, contentMonth),
        () => getNewPipelines(pipelines, contentYear, contentMonth),
        () => getEventsForMonth(events, contentYear, contentMonth),
        () => getEventsForMonth(events, year, month),
        () => getFutureEvents(events, year, month, 1),
        () => getProposalsForMonth(proposals, contentYear, contentMonth),
    ];

    return checks.some((check) => check().length > 0);
}

// The blog / events / advisories collections all share the same "filter to a
// month, sorted by date" and "collect across N adjacent months" shapes. These
// two generics implement that once; each collection only supplies its date
// accessor and sort direction below.
type SortDir = "asc" | "desc";

/** Items whose date falls in the given calendar month, sorted by that date. */
function itemsInMonth<T>(items: T[], getDate: (item: T) => Date, year: number, month: number, dir: SortDir): T[] {
    const sign = dir === "asc" ? 1 : -1;
    return items
        .filter((item) => isInMonth(getDate(item), year, month))
        .sort((a, b) => sign * (getDate(a).getTime() - getDate(b).getTime()));
}

/** Items collected across `count` months adjacent to (year, month), re-sorted as one list. */
function itemsInAdjacentMonths<T>(
    items: T[],
    getDate: (item: T) => Date,
    year: number,
    month: number,
    count: number,
    step: 1 | -1,
    dir: SortDir,
): T[] {
    const sign = dir === "asc" ? 1 : -1;
    return adjacentMonths(year, month, count, step)
        .flatMap(({ year: y, month: m }) => itemsInMonth(items, getDate, y, m, dir))
        .sort((a, b) => sign * (getDate(a).getTime() - getDate(b).getTime()));
}

const blogDate = (post: { data: { pubDate: Date } }) => new Date(post.data.pubDate);
const eventDate = (event: { data: { start: Date } }) => new Date(event.data.start);
const advisoryDate = (advisory: { data: { publishedDate: Date } }) => new Date(advisory.data.publishedDate);

/** Blog posts published in the given month (newest first). */
export function getBlogPostsForMonth<T extends { id: string; data: { pubDate: Date } }>(
    posts: T[],
    year: number,
    month: number,
) {
    return itemsInMonth<T>(posts, blogDate, year, month, "desc");
}

/** Blog posts from the previous N months (for "in case you missed it"). */
export function getBlogPostsForPreviousMonths<T extends { id: string; data: { pubDate: Date } }>(
    posts: T[],
    year: number,
    month: number,
    count: number = 2,
) {
    return itemsInAdjacentMonths<T>(posts, blogDate, year, month, count, -1, "desc");
}

/**
 * Pipeline releases published in the given month.
 * Flags first-ever releases (earliest non-dev release for that pipeline).
 */
export function getPipelineReleasesForMonth(
    pipelines: PipelineWorkflow[],
    year: number,
    month: number,
): NewsletterRelease[] {
    const grouped = new Map<string, NewsletterRelease>();

    for (const pipeline of pipelines) {
        if (pipeline.archived) continue;

        // Find the earliest non-dev release for first-release detection
        const nonDevReleases = pipeline.releases
            .filter((r) => r.tag_name !== "dev")
            .sort((a, b) => new Date(a.published_at).getTime() - new Date(b.published_at).getTime());

        const firstReleaseTag = nonDevReleases.length > 0 ? nonDevReleases[0].tag_name : null;

        for (const release of pipeline.releases) {
            if (release.tag_name === "dev") continue;
            const pubDate = new Date(release.published_at);
            if (!isInMonth(pubDate, year, month)) continue;

            const existing = grouped.get(pipeline.name);
            if (existing) {
                existing.tagNames.push(release.tag_name);
                // Keep the most recent publish date
                if (new Date(release.published_at) > new Date(existing.publishedAt)) {
                    existing.publishedAt = release.published_at;
                }
                // If any release is a first release, flag it
                if (release.tag_name === firstReleaseTag) {
                    existing.isFirstRelease = true;
                }
            } else {
                grouped.set(pipeline.name, {
                    pipelineName: pipeline.name,
                    description: pipeline.description,
                    tagNames: [release.tag_name],
                    publishedAt: release.published_at,
                    isFirstRelease: release.tag_name === firstReleaseTag,
                });
            }
        }
    }

    const releases = [...grouped.values()];

    // Sort tag names within each release by semver ascending (oldest first, left to right)
    for (const release of releases) {
        release.tagNames.sort((a, b) => {
            const aNum = a.replace(/^v/, "").split(".").map(Number);
            const bNum = b.replace(/^v/, "").split(".").map(Number);
            for (let i = 0; i < Math.max(aNum.length, bNum.length); i++) {
                if ((aNum[i] || 0) !== (bNum[i] || 0)) return (aNum[i] || 0) - (bNum[i] || 0);
            }
            return 0;
        });
    }

    // Sort: first releases first, then alphabetically by pipeline name
    releases.sort((a, b) => {
        if (a.isFirstRelease !== b.isFirstRelease) return a.isFirstRelease ? -1 : 1;
        return a.pipelineName.localeCompare(b.pipelineName);
    });

    return releases;
}

/**
 * Pipelines created in the given month (new repos).
 */
export function getNewPipelines(pipelines: PipelineWorkflow[], year: number, month: number): NewPipeline[] {
    return pipelines
        .filter((p) => {
            if (p.archived) return false;
            return isInMonth(new Date(p.created_at), year, month);
        })
        .map((p) => ({
            name: p.name,
            description: p.description,
            createdAt: p.created_at,
            isDevOnly: p.releases.every((r) => r.tag_name === "dev"),
        }))
        .sort((a, b) => a.name.localeCompare(b.name));
}

/** Advisories published in the given month (newest first). */
export function getAdvisoriesForMonth<T extends { id: string; data: { publishedDate: Date } }>(
    advisories: T[],
    year: number,
    month: number,
) {
    return itemsInMonth<T>(advisories, advisoryDate, year, month, "desc");
}

/** Advisories from the previous N months (for "in case you missed it"). */
export function getAdvisoriesForPreviousMonths<T extends { id: string; data: { publishedDate: Date } }>(
    advisories: T[],
    year: number,
    month: number,
    count: number = 2,
) {
    return itemsInAdjacentMonths<T>(advisories, advisoryDate, year, month, count, -1, "desc");
}

/** Events starting in the given month (chronological). */
export function getEventsForMonth<T extends { id: string; data: { start: Date } }>(
    events: T[],
    year: number,
    month: number,
) {
    return itemsInMonth<T>(events, eventDate, year, month, "asc");
}

/** Events starting in the next N months after the given month (chronological). */
export function getFutureEvents<T extends { id: string; data: { start: Date } }>(
    events: T[],
    year: number,
    month: number,
    monthsAhead: number = 2,
) {
    return itemsInAdjacentMonths<T>(events, eventDate, year, month, monthsAhead, 1, "asc");
}

/**
 * Fetch approved proposals from nf-core/proposals GitHub repo.
 * Returns all proposals, to be filtered per-month in the page.
 */
export async function fetchAllProposals(): Promise<RawProposal[]> {
    try {
        const issues = await octokit.paginate(octokit.rest.issues.listForRepo, {
            owner: "nf-core",
            repo: "proposals",
            state: "all",
            per_page: 100,
            headers: {
                "X-GitHub-Api-Version": "2022-11-28",
            },
        });

        return issues
            .filter((issue) => !issue.pull_request) // Exclude PRs
            .map((issue) => ({
                title: issue.title,
                url: issue.html_url,
                number: issue.number,
                labels: issue.labels.map((l) => (typeof l === "string" ? l : l.name || "")),
                closedAt: issue.closed_at,
                stateReason: (issue as any).state_reason ?? null,
                createdAt: issue.created_at || "",
            }));
    } catch (error) {
        console.error("Failed to fetch proposals from nf-core/proposals:", error);
        return [];
    }
}

/**
 * Categorize a proposal title: strip "New pipeline:" prefix for pipeline proposals,
 * keep full title with prefix for others (RFCs, SIGs, etc).
 */
function categorizeProposal(title: string): { category: "pipeline" | "other"; displayTitle: string } {
    const pipelineMatch = title.match(/^New pipeline:\s*(.+)$/i);
    if (pipelineMatch) {
        return { category: "pipeline", displayTitle: pipelineMatch[1].trim() };
    }
    return { category: "other", displayTitle: title };
}

/**
 * Filter proposals to those opened or closed/accepted in a given month.
 * Assigns status ("new" for opened this month, "accepted" for closed this month).
 * A proposal can appear as both if opened and closed in the same month.
 */
export function getProposalsForMonth(proposals: RawProposal[], year: number, month: number): Proposal[] {
    const results = new Map<string, Proposal>();

    for (const p of proposals) {
        const { category, displayTitle } = categorizeProposal(p.title);
        const key = `${p.number}`;

        // Opened this month
        if (isInMonth(new Date(p.createdAt), year, month)) {
            results.set(key, { ...p, displayTitle, category, status: "new" });
        }

        // Closed this month as "completed" = accepted (skip "not_planned" = rejected)
        if (p.closedAt && p.stateReason === "completed" && isInMonth(new Date(p.closedAt), year, month)) {
            const existing = results.get(key);
            if (existing) {
                existing.status = "accepted";
            } else {
                results.set(key, { ...p, displayTitle, category, status: "accepted" });
            }
        }
    }

    return [...results.values()];
}

// ========================================
// SHARED DATA HELPERS
// ========================================

/**
 * Shared getStaticPaths data for newsletter pages: the list of months plus all
 * proposals (fetched once).
 */
export async function getNewsletterStaticPathsData(
    getCollectionFn: (name: string) => Promise<any[]>,
    pipelines: PipelineWorkflow[],
) {
    let blogPosts = await getCollectionFn("blog");
    blogPosts = blogPosts.filter((post: any) => new Date(post.data.pubDate) < new Date());

    let events = await getCollectionFn("events");
    events = events.filter((e: any) => e.id.split("/").length === 2);

    const advisories = await getCollectionFn("advisories");

    let allProposals: RawProposal[] = [];
    try {
        allProposals = await fetchAllProposals();
    } catch (e) {
        console.warn("Could not fetch proposals:", e);
    }

    // Drop any candidate month whose newsletter would render no content (see
    // newsletterMonthHasContent) so we never publish an empty page.
    const months = getNewsletterMonths(blogPosts, events, pipelines, advisories).filter((m) =>
        newsletterMonthHasContent(blogPosts, events, pipelines, advisories, allProposals, m.year, m.month),
    );

    return { months, allProposals };
}

/**
 * The canonical list of newsletter months that actually have content, shared by the
 * listing page and the RSS feed so they stay in lock-step with the generated pages.
 */
export async function getPublishedNewsletterMonths(
    getCollectionFn: (name: string) => Promise<any[]>,
    pipelines: PipelineWorkflow[],
): Promise<NewsletterMonth[]> {
    const { months } = await getNewsletterStaticPathsData(getCollectionFn, pipelines);
    return months;
}

/**
 * The full Astro getStaticPaths array for every newsletter month route
 * (the page, /email, /simple, /markdown). Each entry carries the month plus
 * its neighbours and the month list — the email/simple/markdown routes simply
 * ignore the nav props they don't use.
 */
export async function getNewsletterStaticPaths(
    getCollectionFn: (name: string) => Promise<any[]>,
    pipelines: PipelineWorkflow[],
) {
    const { months, allProposals } = await getNewsletterStaticPathsData(getCollectionFn, pipelines);
    return months.map(({ year, month }, index) => ({
        params: { year: String(year), month: String(month).padStart(2, "0") },
        props: {
            year,
            month,
            allMonths: months,
            allProposals,
            prevMonth: months[index + 1] || null,
            nextMonth: months[index - 1] || null,
        },
    }));
}
