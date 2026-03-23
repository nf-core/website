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
    tagName: string;
    publishedAt: string;
    isFirstRelease: boolean;
}

export interface NewPipeline {
    name: string;
    description: string;
    createdAt: string;
    isDevOnly: boolean;
}

export interface Proposal {
    title: string;
    url: string;
    number: number;
    labels: string[];
    closedAt: string | null;
    createdAt: string;
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
    const releases: NewsletterRelease[] = [];

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

            releases.push({
                pipelineName: pipeline.name,
                description: pipeline.description,
                tagName: release.tag_name,
                publishedAt: release.published_at,
                isFirstRelease: release.tag_name === firstReleaseTag,
            });
        }
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
export async function fetchAllProposals(): Promise<Proposal[]> {
    try {
        const issues = await octokit.paginate("GET /repos/{owner}/{repo}/issues", {
            owner: "nf-core",
            repo: "proposals",
            state: "all",
            per_page: 100,
            headers: {
                "X-GitHub-Api-Version": "2022-11-28",
            },
        });

        return issues
            .filter((issue: any) => !issue.pull_request) // Exclude PRs
            .map((issue: any) => ({
                title: issue.title,
                url: issue.html_url,
                number: issue.number,
                labels: issue.labels.map((l: any) => (typeof l === "string" ? l : l.name)),
                closedAt: issue.closed_at,
                createdAt: issue.created_at,
            }));
    } catch (error) {
        console.error("Failed to fetch proposals from nf-core/proposals:", error);
        return [];
    }
}

/**
 * Filter proposals to those approved/closed in a given month.
 */
export function getProposalsForMonth(proposals: Proposal[], year: number, month: number): Proposal[] {
    return proposals.filter((p) => {
        // Include if closed in this month (approved proposals get closed)
        if (p.closedAt && isInMonth(new Date(p.closedAt), year, month)) return true;
        // Also include if created this month and has an approval-related label
        const approvalLabels = ["approved", "accepted", "pipeline"];
        const hasApprovalLabel = p.labels.some((l) => approvalLabels.includes(l.toLowerCase()));
        if (hasApprovalLabel && isInMonth(new Date(p.createdAt), year, month)) return true;
        return false;
    });
}
