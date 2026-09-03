// Shared helpers for the newsletter renderings. The styled email/web layout
// (NewsletterLayout.astro), the plain-HTML layout (NewsletterSimpleLayout.astro)
// and the markdown renderer below all use these, so the formatting, URLs and
// event-type data stay consistent across every variant.

function stripExt(id: string): string {
    return id.replace(/\.[^/.]+$/, "");
}

/** Normalise an image source: Astro image objects expose `.src`, raw paths are strings. */
export function imgSrc(src: any): string {
    return typeof src === "string" ? src : src?.src;
}

export function blogUrl(baseUrl: string, id: string): string {
    return `${baseUrl}/blog/${stripExt(id)}`;
}
export function eventUrl(baseUrl: string, id: string): string {
    return `${baseUrl}/events/${stripExt(id)}`;
}
export function advisoryUrl(baseUrl: string, id: string): string {
    return `${baseUrl}/advisories/${stripExt(id)}`;
}
export function pipelineUrl(baseUrl: string, name: string): string {
    return `${baseUrl}/${name}`;
}

export function formatDate(date: Date | string): string {
    return new Date(date).toLocaleDateString("en-GB", { day: "numeric", month: "short", year: "numeric" });
}

export function eventDateRange(start: Date, end: Date): string {
    const s = formatDate(start);
    return start.toDateString() !== end.toDateString() ? `${s} – ${formatDate(end)}` : s;
}

// Single source of truth for event types: label + badge colours (the colours are
// only used by the styled layout; the markdown/plain variants use the label).
export const EVENT_TYPES: Record<string, { label: string; color: string; textColor: string }> = {
    bytesize: { label: "Bytesize", color: "#22ae63", textColor: "#fefefe" },
    hackathon: { label: "Hackathon", color: "#0d6efd", textColor: "#fefefe" },
    talk: { label: "Talk", color: "#0dcaf0", textColor: "#212529" },
    training: { label: "Training", color: "#ffc107", textColor: "#212529" },
};

export const eventTypeLabel: Record<string, string> = Object.fromEntries(
    Object.entries(EVENT_TYPES).map(([type, { label }]) => [type, label]),
);

export function formatAdvisoryTypes(type: string | string[]): string[] {
    const types = Array.isArray(type) ? type : [type];
    return types.map((t) =>
        t
            .split("_")
            .map((w) => w.charAt(0).toUpperCase() + w.slice(1))
            .join(" "),
    );
}

export function sortByStatus(a: { status: string }, b: { status: string }): number {
    return a.status === b.status ? 0 : a.status === "accepted" ? -1 : 1;
}

// Escape characters that would break markdown link text.
function mdText(s: string): string {
    return s.replace(/([\[\]])/g, "\\$1");
}

/**
 * Render a newsletter month as plain markdown (for agents / LLMs). No images.
 * `data` is the object returned by getNewsletterContentData.
 */
export function renderNewsletterMarkdown(data: any): string {
    const { year, monthName, baseUrl, newsletterUrl, previewText } = data;
    const lines: string[] = [];
    const proposalStatus = (p: any) => (p.status === "accepted" ? "Accepted" : "New");

    lines.push(`# nf-core/newsletter, ${monthName} ${year}`, "");
    lines.push(previewText, "");
    lines.push(`[View on the nf-core website](${newsletterUrl})`, "");

    if (data.thisMonthBlogPosts.length) {
        lines.push("## Blog posts", "");
        for (const post of data.thisMonthBlogPosts) {
            lines.push(`- [${mdText(post.data.title)}](${blogUrl(baseUrl, post.id)}) — ${post.data.subtitle}`);
            lines.push(`  _${formatDate(post.data.pubDate)} · by ${post.data.authors.join(", ")}_`);
        }
        lines.push("");
    }

    if (data.monthAdvisories.length) {
        lines.push("## Advisories", "");
        for (const a of data.monthAdvisories) {
            const types = formatAdvisoryTypes(a.data.type).join(", ");
            lines.push(
                `- **${String(a.data.severity).toUpperCase()}** [${mdText(a.data.title)}](${advisoryUrl(baseUrl, a.id)})` +
                    `${types ? ` _(${types})_` : ""} — ${a.data.subtitle}`,
            );
        }
        lines.push("");
    }

    if (data.firstReleases.length) {
        lines.push("## New pipeline first releases 🎉", "");
        lines.push("These pipelines just had their very first release.", "");
        for (const r of data.firstReleases) {
            const tags = r.tagNames.map((t: string) => `v${t}`).join(", ");
            lines.push(
                `- [nf-core/${r.pipelineName}](${pipelineUrl(baseUrl, r.pipelineName)}) ${tags} — ${r.description}`,
            );
        }
        lines.push("");
    }

    if (data.otherReleases.length) {
        lines.push("## Pipeline releases", "");
        for (const r of data.otherReleases) {
            const tags = r.tagNames.map((t: string) => `v${t}`).join(", ");
            lines.push(
                `- [nf-core/${r.pipelineName}](${pipelineUrl(baseUrl, r.pipelineName)}) ${tags} _(${formatDate(r.publishedAt)})_`,
            );
        }
        lines.push("");
    }

    if (data.upcomingEvents.length) {
        lines.push("## Upcoming events", "");
        for (const e of data.upcomingEvents) {
            const label = eventTypeLabel[e.data.type] || e.data.type;
            lines.push(
                `- **${label}** [${mdText(e.data.title)}](${eventUrl(baseUrl, e.id)}) — ${eventDateRange(e.data.start, e.data.end)}`,
            );
            if (e.data.subtitle) lines.push(`  ${e.data.subtitle}`);
        }
        lines.push("");
    }

    if (data.recentEvents.length || data.olderBlogPosts.length || data.olderAdvisories.length) {
        lines.push("## In case you missed it", "");
        if (data.recentEvents.length) {
            lines.push("### Recent events", "");
            for (const e of data.recentEvents) {
                const label = eventTypeLabel[e.data.type] || e.data.type;
                lines.push(
                    `- [${mdText(e.data.title)}](${eventUrl(baseUrl, e.id)}) — ${eventDateRange(e.data.start, e.data.end)} (${label})`,
                );
            }
            lines.push("");
        }
        if (data.olderBlogPosts.length) {
            lines.push("### Earlier blog posts", "");
            for (const post of data.olderBlogPosts) {
                lines.push(
                    `- [${mdText(post.data.title)}](${blogUrl(baseUrl, post.id)}) — ${formatDate(post.data.pubDate)}`,
                );
            }
            lines.push("");
        }
        if (data.olderAdvisories.length) {
            lines.push("### Earlier advisories", "");
            for (const a of data.olderAdvisories) {
                lines.push(
                    `- **${String(a.data.severity).toUpperCase()}** [${mdText(a.data.title)}](${advisoryUrl(baseUrl, a.id)}) — ${formatDate(a.data.publishedDate)}`,
                );
            }
            lines.push("");
        }
    }

    if (data.newPipelines.length || data.pipelineProposals.length || data.otherProposals.length) {
        lines.push("## New pipelines & proposals", "");
        if (data.newPipelines.length) {
            lines.push("### New pipeline repositories", "");
            for (const p of data.newPipelines) {
                lines.push(
                    `- [nf-core/${p.name}](${pipelineUrl(baseUrl, p.name)})${p.isDevOnly ? " _(in development)_" : ""} — ${p.description}`,
                );
            }
            lines.push("");
        }
        if (data.pipelineProposals.length) {
            lines.push("### Pipeline proposals", "");
            lines.push("From [nf-core/proposals](https://github.com/nf-core/proposals/issues).", "");
            for (const p of [...data.pipelineProposals].sort(sortByStatus)) {
                lines.push(`- [${mdText(p.displayTitle)}](${p.url}) — ${proposalStatus(p)} (#${p.number})`);
            }
            lines.push("");
        }
        if (data.otherProposals.length) {
            lines.push("### Other proposals", "");
            for (const p of [...data.otherProposals].sort(sortByStatus)) {
                lines.push(`- [${mdText(p.displayTitle)}](${p.url}) — ${proposalStatus(p)} (#${p.number})`);
            }
            lines.push("");
        }
    }

    lines.push("---", "");
    lines.push("This newsletter is brought to you by the [nf-core](https://nf-co.re) community.", "");
    return lines.join("\n");
}
