import { octokit } from "@components/octokit.js";

// ========================================
// TYPE DEFINITIONS
// ========================================

interface PipelineWorkflow {
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
 * Blog posts published in the given month.
 */
export function getBlogPostsForMonth(
    posts: { id: string; data: { pubDate: Date; title: string; subtitle: string; authors: string[] } }[],
    year: number,
    month: number,
) {
    return posts
        .filter((post) => isInMonth(new Date(post.data.pubDate), year, month))
        .sort((a, b) => new Date(b.data.pubDate).getTime() - new Date(a.data.pubDate).getTime());
}

/**
 * Blog posts from the previous N months (for "in case you missed it").
 */
export function getBlogPostsForPreviousMonths(
    posts: { id: string; data: { pubDate: Date; title: string; subtitle: string; authors: string[] } }[],
    year: number,
    month: number,
    count: number = 2,
) {
    const results: typeof posts = [];
    let y = year;
    let m = month;

    for (let i = 0; i < count; i++) {
        m--;
        if (m === 0) {
            m = 12;
            y--;
        }
        results.push(...getBlogPostsForMonth(posts, y, m));
    }

    return results.sort((a, b) => new Date(b.data.pubDate).getTime() - new Date(a.data.pubDate).getTime());
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

/**
 * Events starting in the given month.
 */
export function getEventsForMonth(
    events: { id: string; data: { start: Date; end: Date; title: string; subtitle: string; type: string } }[],
    year: number,
    month: number,
) {
    return events
        .filter((e) => isInMonth(new Date(e.data.start), year, month))
        .sort((a, b) => new Date(a.data.start).getTime() - new Date(b.data.start).getTime());
}

/**
 * Advisories published in the given month.
 */
export function getAdvisoriesForMonth(
    advisories: {
        id: string;
        data: { publishedDate: Date; title: string; subtitle: string; severity: string; type: string[] };
    }[],
    year: number,
    month: number,
) {
    return advisories
        .filter((a) => isInMonth(new Date(a.data.publishedDate), year, month))
        .sort((a, b) => new Date(b.data.publishedDate).getTime() - new Date(a.data.publishedDate).getTime());
}

/**
 * Events starting in the next N months after the given month.
 */
export function getFutureEvents(
    events: { id: string; data: { start: Date; end: Date; title: string; subtitle: string; type: string } }[],
    year: number,
    month: number,
    monthsAhead: number = 2,
) {
    const results: typeof events = [];
    let y = year;
    let m = month;

    for (let i = 0; i < monthsAhead; i++) {
        m++;
        if (m === 13) {
            m = 1;
            y++;
        }
        results.push(...getEventsForMonth(events, y, m));
    }

    return results.sort((a, b) => new Date(a.data.start).getTime() - new Date(b.data.start).getTime());
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
 * Shared getStaticPaths data for newsletter pages.
 * Both [month].astro and [month]/email.astro call this.
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
    const months = getNewsletterMonths(blogPosts, events, pipelines, advisories);

    let allProposals: any[] = [];
    try {
        allProposals = await fetchAllProposals();
    } catch (e) {
        console.warn("Could not fetch proposals:", e);
    }

    return { months, allProposals };
}

/**
 * Shared content data fetching for a single newsletter month.
 * Returns the data object expected by NewsletterContent component.
 */
export async function getNewsletterContentData(
    getCollectionFn: (name: string) => Promise<any[]>,
    pipelines: PipelineWorkflow[],
    year: number,
    month: number,
    allProposals: any[],
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
    if (releases.length > 0)
        previewParts.push(`${releases.length} pipeline release${releases.length > 1 ? "s" : ""}`);
    if (upcomingEvents.length > 0)
        previewParts.push(`${upcomingEvents.length} upcoming event${upcomingEvents.length > 1 ? "s" : ""}`);
    const previewText =
        previewParts.length > 0 ? previewParts.join(", ") + "." : `Community news for ${monthName} ${year}.`;
    const newsletterUrl = `${baseUrl}/newsletter/${year}/${String(month).padStart(2, "0")}`;

    // Pre-resolve all image sources
    async function resolveImageSrc(headerImage: string): Promise<any> {
        if (headerImage.startsWith("/assets/")) {
            const key = "/src" + headerImage;
            const loader = images[key];
            if (loader) {
                const mod = (await loader()) as any;
                return mod.default;
            }
        }
        return headerImage;
    }

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
        blogImageSrcs,
        eventImageSrcs,
    };
}
