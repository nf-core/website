// GitHub-markdown transformations for pipeline READMEs/docs, rendered via
// renderPipelineMarkdown (sites/main-site/utils/pipelineMarkdown.ts). Originally
// a remark plugin (bin/remark-github-markdown.js, since removed); this is the
// satteri port that replaced it.
// - Fixes image URLs
// - Fixes link URLs
// - Rewrites theme-aware <picture> into hide-dark/hide-light <img> pairs
// - Converts GitHub admonitions ([!NOTE] etc.) into directives for admonitions.ts
// - Cleans up headers
// - Fixes code blocks
//
// The README trimming (cut everything before "## Introduction", drop ":warning:"
// sections) worked on root.children in the old remark plugin; satteri has no root
// visitor, so it is done on the markdown source instead — see trimGitHubReadme below.
import { defineMdastPlugin } from "satteri";

const SPECIAL_FILES_REGEX = /^(\.github\/CONTRIBUTING\.md|CITATIONS\.md|CHANGELOG\.md|\.config)$/;
const MDX_ANCHOR_REGEX = /\.mdx?#/;
const ADMONITION_REGEX = /^\[!(NOTE|TIP|IMPORTANT|WARNING|CAUTION)\]\s*(.*)/;
const IMG_SRC_REGEX = /<img(.*?)src="(.*?)"/g;
const SOURCE_SRCSET_REGEX = /<source(.*?)srcset="(.*?)"/g;
const PICTURE_REGEX = /<picture\b[^>]*>([\s\S]*?)<\/picture>/gi;
const SOURCE_TAG_REGEX = /<source\b([^>]*?)\/?>/gi;
const IMG_TAG_REGEX = /<img\b([^>]*?)\/?>/i;

function getHtmlAttr(attrs: string, name: string): string | undefined {
    const match = attrs.match(new RegExp(`\\b${name}\\s*=\\s*"([^"]*)"`, "i"));
    return match ? match[1] : undefined;
}

/**
 * Rewrite GitHub-style theme-aware <picture> elements into the dual-<img>
 * `hide-dark`/`hide-light` pattern used elsewhere on the website (see
 * DarkModeImage.astro), so they follow the site theme switcher (data-bs-theme)
 * via CSS rather than only the OS `prefers-color-scheme`.
 *
 * A native <picture> picks its <source> from the *OS* colour scheme, which the
 * website's theme switcher can't override without client JS. On GitHub the
 * plugin doesn't run, so the original <picture> is left intact and still works.
 *
 * Only <picture> elements whose <source> uses a prefers-color-scheme media
 * query are converted; width-based art direction is left untouched.
 */
function rewriteThemedPictures(value: string): string {
    return value.replace(PICTURE_REGEX, (whole: string, inner: string) => {
        // Collect srcsets keyed by the prefers-color-scheme they target.
        const sources: Record<string, string> = {};
        for (const [, attrs] of inner.matchAll(SOURCE_TAG_REGEX)) {
            const srcset = getHtmlAttr(attrs, "srcset");
            const scheme = (getHtmlAttr(attrs, "media") || "").match(/prefers-color-scheme\s*:\s*(dark|light)/i)?.[1];
            if (srcset && scheme) sources[scheme.toLowerCase()] = srcset;
        }
        const imgAttrs = inner.match(IMG_TAG_REGEX)?.[1];
        if (imgAttrs === undefined || (!sources.dark && !sources.light)) return whole;

        // The fallback <img> fills whichever variant has no explicit <source>.
        const fallback = getHtmlAttr(imgAttrs, "src");
        const light = sources.light ?? fallback;
        const dark = sources.dark ?? fallback;
        if (!light || !dark) return whole;

        // Carry over the fallback <img>'s other attributes (alt, width, …),
        // setting src/class explicitly per variant.
        const cls = getHtmlAttr(imgAttrs, "class");
        const base = imgAttrs.replace(/\b(?:src|class)\s*=\s*"[^"]*"/gi, "").replace(/\s+/g, " ").trim();
        const img = (src: string, hide: string) =>
            `<img ${base ? base + " " : ""}src="${src}" class="${[hide, cls].filter(Boolean).join(" ")}" />`;

        return img(light, "hide-dark") + img(dark, "hide-light");
    });
}

export interface GitHubMarkdownOptions {
    org?: string;
    repo: string;
    ref: string;
    parentDirectory?: string;
}

export function createGitHubMarkdownPlugin(options: GitHubMarkdownOptions) {
    const { org = "nf-core", repo, ref, parentDirectory = "" } = options;
    const baseRawUrl = `https://raw.githubusercontent.com/${org}/${repo}/${ref}/`;
    const baseRepoUrl = `https://github.com/${org}/${repo}/blob/${ref}/`;

    return defineMdastPlugin({
        name: "github-markdown",
        image(node: any, ctx) {
            if (node.url && !node.url.startsWith("http")) {
                ctx.setProperty(node, "url", `${baseRawUrl}${parentDirectory}/${node.url}`);
            }
        },
        link(node: any, ctx) {
            if (!node.url || node.url.startsWith("http")) return;
            if (SPECIAL_FILES_REGEX.test(node.url)) {
                ctx.setProperty(node, "url", `${baseRepoUrl}${node.url}`);
            } else if (node.url.includes("../assets/")) {
                ctx.setProperty(node, "url", `${baseRepoUrl}${node.url.replace("../assets/", "assets/")}`);
            } else if (MDX_ANCHOR_REGEX.test(node.url)) {
                ctx.setProperty(node, "url", node.url.replace(MDX_ANCHOR_REGEX, "#"));
            }
        },
        code(node: any, ctx) {
            if (node.lang === "nextflow") {
                ctx.setProperty(node, "lang", "groovy");
            }
        },
        heading(node: any, ctx) {
            if (node.depth !== 1 || !node.children) return;
            const headingText = ctx.textContent(node);
            const prefix = `nf-core/${repo}: `;
            if (headingText.startsWith(prefix) && node.children[0]?.type === "text") {
                const children = [...node.children];
                children[0] = { ...children[0], value: children[0].value.replace(prefix, "") };
                return { ...node, children };
            }
        },
        blockquote(node: any) {
            // GitHub-style admonitions: > [!NOTE] ... → :::note directive, which
            // admonitions.ts (running after this plugin) turns into the styled box.
            const firstChild = node.children?.[0];
            if (firstChild?.type !== "paragraph" || !firstChild.children?.length) return;
            const firstText = firstChild.children[0];
            if (firstText.type !== "text" || !firstText.value.startsWith("[!")) return;
            const match = firstText.value.match(ADMONITION_REGEX);
            if (!match) return;

            let type = match[1].toLowerCase();
            if (type === "important") type = "info";
            else if (type === "caution") type = "danger";

            const attributes: Record<string, string> = {};
            if (match[1] === "IMPORTANT") attributes.title = "Important";
            else if (match[1] === "CAUTION") attributes.title = "Caution";

            const paragraphChildren = [
                ...(match[2] ? [{ type: "text", value: match[2] + " " }] : []),
                ...firstChild.children.slice(1),
            ];

            return {
                type: "containerDirective",
                name: type,
                attributes,
                children: [
                    ...(paragraphChildren.length > 0 ? [{ type: "paragraph", children: paragraphChildren }] : []),
                    ...node.children.slice(1),
                ],
            };
        },
        html(node: any, ctx) {
            let value = node.value;
            if (value.includes("<img") && !value.includes('src="http')) {
                value = value.replace(
                    IMG_SRC_REGEX,
                    (_match: string, attrs: string, src: string) =>
                        `<img${attrs}src="${baseRawUrl}${parentDirectory}/${src}"`,
                );
            }
            if (value.includes("<source") && !value.includes('srcset="http')) {
                value = value.replace(
                    SOURCE_SRCSET_REGEX,
                    (_match: string, attrs: string, src: string) =>
                        `<source${attrs}srcset="${baseRawUrl}${parentDirectory}/${src}"`,
                );
            }
            // Convert theme-aware <picture> elements (now with absolute srcset/src)
            // into the CSS hide-dark/hide-light pattern so they follow the site theme.
            if (value.includes("<picture")) {
                value = rewriteThemedPictures(value);
            }
            if (value !== node.value) {
                ctx.setProperty(node, "value", value);
            }
        },
    });
}

interface HeadingLine {
    line: number;
    depth: number;
    text: string;
}

function findHeadings(lines: string[]): HeadingLine[] {
    const headings: HeadingLine[] = [];
    let inFence = false;
    lines.forEach((line, i) => {
        if (/^\s*(```|~~~)/.test(line)) {
            inFence = !inFence;
            return;
        }
        if (inFence) return;
        const match = line.match(/^(#{1,6})\s+(.*)/);
        if (match) headings.push({ line: i, depth: match[1].length, text: match[2].trim() });
    });
    return headings;
}

/**
 * README cleanup, previously the tree-level passes of remark-github-markdown:
 * drop everything before the first "## Introduction" heading and remove
 * ":warning:" sections (warning heading up to the next heading).
 */
export function trimGitHubReadme(source: string): string {
    let lines = source.split("\n");
    let headings = findHeadings(lines);

    const intro = headings.find((h) => h.depth === 2 && h.text === "Introduction");
    if (intro && intro.line > 0) {
        lines = lines.slice(intro.line);
        headings = headings.filter((h) => h.line >= intro.line).map((h) => ({ ...h, line: h.line - intro.line }));
    }

    const warning = headings.find((h) => h.depth === 2 && h.text.includes(":warning:"));
    if (warning) {
        const next = headings.find((h) => h.line > warning.line);
        if (next) {
            lines.splice(warning.line, next.line - warning.line);
        }
    }

    return lines.join("\n");
}
