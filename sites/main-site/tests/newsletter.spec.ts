import { test, expect } from "@playwright/test";
import {
    getMonthName,
    getNewsletterMonths,
    getBlogPostsForMonth,
    getBlogPostsForPreviousMonths,
    getPipelineReleasesForMonth,
    getProposalsForMonth,
    newsletterMonthHasContent,
} from "@utils/newsletter";

// Unit tests for the newsletter data/date logic. These functions are pure and
// build-time critical: the month maths (incl. year boundaries), release grouping
// and proposal categorisation are easy to break silently and invisible in review.

const blog = (id: string, pubDate: string) => ({ id, data: { pubDate: new Date(pubDate) } });
const advisory = (id: string, publishedDate: string) => ({ id, data: { publishedDate: new Date(publishedDate) } });
const evt = (id: string, start: string) => ({ id, data: { start: new Date(start) } });

test("getMonthName maps 1-indexed months to English names", () => {
    expect(getMonthName(1)).toBe("January");
    expect(getMonthName(12)).toBe("December");
});

test("getNewsletterMonths dedupes, sorts newest-first, and drops future months", () => {
    const future = new Date();
    future.setFullYear(future.getFullYear() + 5);

    const months = getNewsletterMonths(
        [blog("a", "2020-03-10"), blog("future", future.toISOString())],
        [evt("e", "2020-05-20")],
        [], // pipelines
        [advisory("adv", "2020-05-02")], // same month as the event -> must dedupe
    );

    // 2020-05 appears once despite event + advisory; future month excluded.
    expect(months).toContainEqual({ year: 2020, month: 5 });
    expect(months).toContainEqual({ year: 2020, month: 3 });
    expect(months.filter((m) => m.year === 2020 && m.month === 5)).toHaveLength(1);
    expect(months.some((m) => m.year === future.getFullYear())).toBe(false);

    // Sorted strictly newest-first.
    const ordinals = months.map((m) => m.year * 12 + m.month);
    expect(ordinals).toEqual([...ordinals].sort((a, b) => b - a));
});

test("newsletterMonthHasContent gates on the content the newsletter actually renders", () => {
    // A newsletter dated the 1st of a month looks back at the *previous* calendar
    // month. A blog post in March 2020 therefore belongs in the April 2020 newsletter,
    // not March's.
    const blogPosts = [blog("a", "2020-03-10")];

    // April 2020 shows March's content -> has content.
    expect(newsletterMonthHasContent(blogPosts, [], [], [], [], 2020, 4)).toBe(true);
    // March 2020 shows February's content (empty) -> no content.
    expect(newsletterMonthHasContent(blogPosts, [], [], [], [], 2020, 3)).toBe(false);

    // A lone advisory dated May 2020 surfaces in the June 2020 newsletter, leaving
    // May's own page empty (the bug behind the orphaned /newsletter/2001/05 page).
    const advisories = [advisory("adv", "2020-05-25")];
    expect(newsletterMonthHasContent([], [], [], advisories, [], 2020, 5)).toBe(false);
    expect(newsletterMonthHasContent([], [], [], advisories, [], 2020, 6)).toBe(true);

    // Upcoming events count for the current and next month's newsletter.
    const events = [evt("e", "2020-05-20")];
    expect(newsletterMonthHasContent([], events, [], [], [], 2020, 5)).toBe(true); // event this month
    expect(newsletterMonthHasContent([], events, [], [], [], 2020, 4)).toBe(true); // event next month
});

test("getBlogPostsForMonth keeps only that month, newest first", () => {
    const posts = [blog("old", "2025-03-02"), blog("new", "2025-03-28"), blog("other", "2025-04-01")];
    const result = getBlogPostsForMonth(posts, 2025, 3);
    expect(result.map((p) => p.id)).toEqual(["new", "old"]);
});

test("getBlogPostsForPreviousMonths walks back across the year boundary", () => {
    // For January 2025 the previous two months are Dec 2024 and Nov 2024.
    const posts = [blog("dec", "2024-12-15"), blog("nov", "2024-11-10"), blog("jan", "2025-01-05")];
    const result = getBlogPostsForPreviousMonths(posts, 2025, 1, 2);
    const ids = result.map((p) => p.id);
    expect(ids).toContain("dec");
    expect(ids).toContain("nov");
    expect(ids).not.toContain("jan"); // the current month is not "previous"
});

test("getPipelineReleasesForMonth groups tags, flags first releases, and sorts tags by semver", () => {
    const pipelines = [
        {
            name: "demo",
            description: "Demo pipeline",
            created_at: "2024-01-01T00:00:00Z",
            archived: false,
            topics: [],
            stargazers_count: 0,
            releases: [
                { tag_name: "dev", published_at: "2025-03-10", tag_sha: "", has_schema: false },
                { tag_name: "1.10.0", published_at: "2025-03-20", tag_sha: "", has_schema: false },
                { tag_name: "1.9.0", published_at: "2025-03-01", tag_sha: "", has_schema: false },
            ],
        },
    ];
    const releases = getPipelineReleasesForMonth(pipelines as any, 2025, 3);

    expect(releases).toHaveLength(1);
    expect(releases[0].pipelineName).toBe("demo");
    // "dev" excluded; tags sorted numerically (1.9.0 before 1.10.0, not lexically).
    expect(releases[0].tagNames).toEqual(["1.9.0", "1.10.0"]);
    // 1.9.0 is the earliest non-dev release for this pipeline, so this is a first release.
    expect(releases[0].isFirstRelease).toBe(true);
    // publishedAt tracks the most recent release in the group.
    expect(releases[0].publishedAt).toBe("2025-03-20");
});

test("getProposalsForMonth categorises and assigns status", () => {
    const proposals = [
        // Opened this month -> new pipeline proposal (prefix stripped).
        {
            title: "New pipeline: foobar",
            url: "u1",
            number: 1,
            labels: [],
            closedAt: null,
            stateReason: null,
            createdAt: "2025-06-10",
        },
        // Closed as completed this month (opened earlier) -> accepted, non-pipeline.
        {
            title: "RFC: governance",
            url: "u2",
            number: 2,
            labels: [],
            closedAt: "2025-06-20",
            stateReason: "completed",
            createdAt: "2025-01-01",
        },
        // Closed as not_planned -> not accepted; still "new" because opened this month.
        {
            title: "New pipeline: rejected",
            url: "u3",
            number: 3,
            labels: [],
            closedAt: "2025-06-21",
            stateReason: "not_planned",
            createdAt: "2025-06-05",
        },
    ];
    const result = getProposalsForMonth(proposals, 2025, 6);
    const byNumber = Object.fromEntries(result.map((p) => [p.number, p]));

    expect(byNumber[1]).toMatchObject({ category: "pipeline", displayTitle: "foobar", status: "new" });
    expect(byNumber[2]).toMatchObject({ category: "other", displayTitle: "RFC: governance", status: "accepted" });
    expect(byNumber[3]).toMatchObject({ category: "pipeline", displayTitle: "rejected", status: "new" });
});
