// The docs of the nf-core-utils plugin for /nf-core-utils/: fetch them from GitHub, join them on one page and adapt
// their MkDocs and GitHub markdown to the website.
import GithubSlugger from "github-slugger";

export const docsRawUrl = "https://raw.githubusercontent.com/nf-core/nf-core-utils/refs/heads/main/docs/";
export const docsBlobUrl = "https://github.com/nf-core/nf-core-utils/blob/main/docs/";
export const docsEntry = "NfCoreUtilities.md";

export interface NfCoreUtilsDoc {
    file: string; // path relative to docs/
    content: string;
}

export interface TocHeading {
    text: string; // HTML
    slug: string;
    depth: number;
}

async function fetchDoc(file: string): Promise<string | null> {
    try {
        const res = await fetch(docsRawUrl + file);
        if (!res.ok) {
            throw new Error(`HTTP ${res.status}`);
        }
        return (await res.text()).replace(/\r\n?/g, "\n");
    } catch (error) {
        console.warn(`nf-core-utils docs: could not fetch ${file}: ${error}`);
        return null;
    }
}

// Start from the guide and follow its links to the docs of the extensions and utilities, the other files in docs/ are
// for the developers of the plugin. Docs at the top of docs/ come first, then those of each subfolder.
export async function fetchNfCoreUtilsDocs(): Promise<NfCoreUtilsDoc[]> {
    const docs: NfCoreUtilsDoc[] = [];
    const seen = new Set([docsEntry]);
    let files = [docsEntry];
    while (files.length) {
        const contents = await Promise.all(files.map(fetchDoc));
        const next: string[] = [];
        files.forEach((file, index) => {
            const content = contents[index];
            // a link to a doc that could not be fetched keeps pointing to GitHub
            if (content === null) {
                return;
            }
            docs.push({ file, content });
            for (const target of linkedDocs(file, content)) {
                if (!seen.has(target)) {
                    seen.add(target);
                    next.push(target);
                }
            }
        });
        files = next;
    }
    return docs.sort((a, b) => Number(a.file.includes("/")) - Number(b.file.includes("/")));
}

function linkedDocs(file: string, content: string): string[] {
    return splitCode(content).flatMap((item) =>
        item.code
            ? []
            : links(item.line)
                  .filter((link) => !link.image)
                  .map((link) => new URL(link.target, docsRawUrl + file))
                  .filter((url) => url.href.startsWith(docsRawUrl) && url.pathname.endsWith(".md"))
                  .map((url) => decodeURIComponent(url.pathname.slice(new URL(docsRawUrl).pathname.length))),
    );
}

// --- Code blocks

type Item = { code: false; line: string } | { code: true; lines: string[] };

// Split markdown into fenced code blocks and the other lines. A fence closes with a fence of the same character, at
// least as long and without info string. A fence still open at the end is closed, not to swallow what comes after.
export function splitCode(md: string): Item[] {
    const items: Item[] = [];
    let block: string[] | undefined;
    let fence = "";
    for (const line of md.split("\n")) {
        if (block) {
            block.push(line);
            const close = line.match(/^ *(`{3,}|~{3,})[ \t]*$/);
            if (close && close[1][0] === fence[0] && close[1].length >= fence.length) {
                items.push({ code: true, lines: block });
                block = undefined;
            }
            continue;
        }
        const open = line.match(/^( *)(`{3,}|~{3,})(.*)$/);
        // the info string of a backtick fence has no backtick, else it is inline code
        if (open && !(open[2][0] === "`" && open[3].includes("`"))) {
            fence = open[2];
            block = [line];
            continue;
        }
        items.push({ code: false, line });
    }
    if (block) {
        block.push(block[0].match(/^ */)![0] + fence);
        items.push({ code: true, lines: block });
    }
    return items;
}

function joinItems(items: Item[]): string {
    return items.flatMap((item) => (item.code ? item.lines : [item.line])).join("\n");
}

// --- Inline code and links

// [start, end] of the code spans of a line
function codeSpans(line: string): [number, number][] {
    return [...line.matchAll(/(`+)(?!`)[\s\S]*?[^`]\1(?!`)|(`+)(?!`)\S\2(?!`)/g)].map((match) => [
        match.index!,
        match.index! + match[0].length,
    ]);
}

interface Link {
    start: number;
    end: number;
    image: boolean;
    text: string;
    target: string;
    title: string;
}

// inline links and images outside of code spans, `[text](target "title")`
function links(line: string): Link[] {
    const spans = codeSpans(line);
    return [...line.matchAll(/(!?)\[((?:[^\][]|\[[^\]]*\])*)\]\(\s*<?([^)\s>]*)>?((?:\s+(?:"[^"]*"|'[^']*'))?)\s*\)/g)]
        .filter((match) => !spans.some(([start, end]) => match.index! > start && match.index! < end))
        .map((match) => ({
            start: match.index!,
            end: match.index! + match[0].length,
            image: match[1] === "!",
            text: match[2],
            target: match[3],
            title: match[4],
        }));
}

function replaceLinks(line: string, replace: (link: Link) => string): string {
    let out = "";
    let last = 0;
    for (const link of links(line)) {
        out += line.slice(last, link.start) + replace(link);
        last = link.end;
    }
    return out + line.slice(last);
}

// --- Emojis

const variationSelector = String.fromCodePoint(0xfe0f);
const zeroWidthJoiner = String.fromCodePoint(0x200d);
// characters shown as emojis, not symbols like ↔, ℹ, ©, ® and ™ that are only emojis with a variation selector
const emojiCharacter = `(?:\\p{Emoji_Presentation}|\\p{Extended_Pictographic}${variationSelector})\\p{Emoji_Modifier}?`;
const emoji = `${emojiCharacter}(?:${zeroWidthJoiner}(?:${emojiCharacter}|\\p{Extended_Pictographic}))*${variationSelector}? ?`;

// emojis, with the space after them, outside of code spans
function stripEmojis(line: string): string {
    return line.replace(new RegExp(`(\`+)(?!\`)[\\s\\S]*?[^\`]\\1(?!\`)|${emoji}`, "gu"), (match, ticks) =>
        ticks ? match : "",
    );
}

// emojis marking code comments, as in `// ✅ CORRECT`, the ones in strings can be part of the example
function stripCommentEmojis(line: string): string {
    return line.replace(new RegExp(`(^|\\s)(\\/\\/|#)(\\s*)${emoji}`, "gu"), "$1$2$3");
}

// --- Admonitions

// MkDocs and GitHub admonition types to the website ones
const admonitionTypes: Record<string, string> = {
    note: "note",
    abstract: "note",
    info: "info",
    important: "info",
    question: "info",
    example: "info",
    tip: "tip",
    success: "tip",
    warning: "warning",
    caution: "danger",
    danger: "danger",
    failure: "danger",
    bug: "danger",
};

function admonition(indent: string, type: string, title: string | undefined, body: string[]): string[] {
    // titles are text in a directive attribute
    const clean = stripEmojis(title ?? "")
        .replace(/`/g, "")
        .replace(/"/g, "'")
        .replace(/\{/g, "(")
        .replace(/\}/g, ")")
        .trim();
    while (body.length && !body[0].trim()) {
        body.shift();
    }
    while (body.length && !body[body.length - 1].trim()) {
        body.pop();
    }
    return [
        `${indent}:::${admonitionTypes[type.toLowerCase()] ?? "note"}${clean ? `{title="${clean}"}` : ""}`,
        ...toWebsiteMarkdown(body.join("\n"))
            .split("\n")
            .map((line) => (line ? indent + line : line)),
        `${indent}:::`,
    ];
}

const dedent = (line: string, width: number) => line.replace(new RegExp(`^ {0,${width}}`), "");

// Outside of code blocks:
// - turn MkDocs admonitions into website ones (`!!! tip "Title"` → `:::tip{title="Title"}`), with as content either
//   the indented blocks after it, the paragraph right after it, or a `- a - b` list on the same line
// - turn GitHub alerts (`> [!TIP]`) into website admonitions
// - remove emojis
// and in code blocks, remove the emojis marking comments
export function toWebsiteMarkdown(md: string): string {
    const items = splitCode(md);
    const out: string[] = [];
    for (let i = 0; i < items.length; i++) {
        const item = items[i];
        if (item.code) {
            out.push(...item.lines.map(stripCommentEmojis));
            continue;
        }
        const alert = item.line.match(/^( *)> *\[!(\w+)\][ \t]*(.*)$/);
        if (alert) {
            const [, indent, type, title] = alert;
            const quote = new RegExp(`^${indent}> ?`);
            const body: string[] = [];
            // the quoted lines, and the lines continuing their paragraph
            while (i + 1 < items.length && !items[i + 1].code) {
                const line = (items[i + 1] as { line: string }).line;
                const lazy =
                    body.length &&
                    body[body.length - 1].trim() &&
                    line.trim() &&
                    !/^ *([#>|]|[-*+] |\d+\. )/.test(line);
                if (!quote.test(line) && !lazy) {
                    break;
                }
                body.push(line.replace(quote, ""));
                i++;
            }
            out.push(...admonition(indent, type, title.replace(/^["“'](.*)["”']$/, "$1"), body));
            continue;
        }
        const mkdocs = item.line.match(
            /^( *)(?:!!!|\?\?\?\+?)[ \t]*(\w+)(?:[ \t]+(?:"([^"]*)"|'([^']*)'|“([^”]*)”))?(.*)$/,
        );
        if (!mkdocs) {
            out.push(stripEmojis(item.line));
            continue;
        }
        const [, indent, type, double, single, curly, rest] = mkdocs;
        const inline = rest.trim();
        const body = inline.startsWith("- ")
            ? inline
                  .slice(2)
                  .split(/ - (?=\S)/)
                  .map((entry) => `- ${entry.trim()}`)
            : inline
              ? [inline]
              : [];
        // MkDocs content: the blocks indented by 4 more spaces, blank lines included
        const width = indent.length + 4;
        const indented = (next: Item) =>
            (next.code ? next.lines[0] : next.line).startsWith(" ".repeat(width)) || (!next.code && !next.line.trim());
        let next = i + 1;
        while (next < items.length && !items[next].code && !(items[next] as { line: string }).line.trim()) {
            next++;
        }
        if (
            next < items.length &&
            (items[next].code ? items[next].lines[0] : items[next].line).startsWith(" ".repeat(width))
        ) {
            while (i + 1 < items.length && indented(items[i + 1])) {
                const block = items[++i];
                body.push(...(block.code ? block.lines : [block.line]).map((line) => dedent(line, width)));
            }
            // the blank lines after the content stay after the admonition
            while (!items[i].code && !(items[i] as { line: string }).line.trim()) {
                body.pop();
                i--;
            }
        } else if (!body.length) {
            // the paragraph right after it, with its code blocks
            while (i + 1 < items.length && (items[i + 1].code || (items[i + 1] as { line: string }).line.trim())) {
                const block = items[++i];
                body.push(...(block.code ? block.lines : [block.line]).map((line) => dedent(line, indent.length)));
            }
        }
        out.push(...admonition(indent, type, double ?? single ?? curly, body));
    }
    return out.join("\n");
}

// --- Headings

interface Heading {
    item: number; // index of the line in the items
    depth: number;
    text: string;
}

// ATX and setext headings outside of code blocks
function headings(items: Item[]): Heading[] {
    const found: Heading[] = [];
    items.forEach((item, index) => {
        if (item.code) {
            return;
        }
        const atx = item.line.match(/^ {0,3}(#{1,6})(?:[ \t]+(.*?))??(?:[ \t]+#+)?[ \t]*$/);
        if (atx) {
            found.push({ item: index, depth: atx[1].length, text: atx[2] ?? "" });
            return;
        }
        const next = items[index + 1];
        const previous = items[index - 1];
        const underline = next && !next.code && next.line.match(/^ {0,3}(=+|-+)[ \t]*$/);
        const paragraph =
            item.line.trim() &&
            !/^ {0,3}([>|<]|[-*+] |\d+[.)] |:::|!!!)/.test(item.line) &&
            (!previous || previous.code || !previous.line.trim());
        if (underline && paragraph) {
            found.push({ item: index, depth: underline[1][0] === "=" ? 1 : 2, text: item.line.trim() });
        }
    });
    return found;
}

// the heading text without its markdown, as GitHub makes its anchor from the rendered text
function plainText(text: string): string {
    return text
        .replace(/!?\[((?:[^\][]|\[[^\]]*\])*)\]\([^)]*\)/g, "$1")
        .replace(/<[^>]+>/g, "")
        .replace(/`/g, "")
        .replace(/(^|[\s(])[*_]+|[*_]+(?=[\s).,:;!?]|$)/g, "$1");
}

// --- The page

const marker = (index: number) => `<!-- nf-core-utils doc ${index} -->`;
const placeholder = (index: number, anchor: string) => `#nf-core-utils-${index}-${anchor}`;

// The docs joined in one markdown, for the website. The docs number their sections themselves, from 1 in each doc and
// not all of them: number them across the page instead, each doc as a chapter (`4.`) and its sections (`4.1.`).
// Links to the docs are placeholders resolved in the rendered HTML by resolveNfCoreUtilsHtml, as the ids of the
// headings are only known once rendered.
export function joinNfCoreUtilsDocs(docs: NfCoreUtilsDoc[]): string {
    if (!docs.length) {
        return "# Error loading content\n\nUnable to load the documentation from GitHub.";
    }
    let chapter = 0;
    return docs
        .map((doc, index) => {
            const items = splitCode(toWebsiteMarkdown(doc.content));
            // the section numbers of the doc, to link "section 4.2" to its heading
            const sections = new Map<string, { heading: number; title: string }>();
            let section = 0;
            headings(items).forEach((heading, number) => {
                const [, old, title] = heading.text.match(/^(?:(\d+(?:\.\d+)*)\.[ \t]+)?(.*)$/)!;
                if (old) {
                    sections.set(old, { heading: number, title: plainText(title) });
                }
                if (heading.depth === 1) {
                    chapter++;
                    section = 0;
                }
                const prefix =
                    heading.depth === 1 ? `${chapter}. ` : heading.depth === 2 ? `${chapter}.${++section}. ` : "";
                items[heading.item] = { code: false, line: `${"#".repeat(heading.depth)} ${prefix}${title}` };
                // the underline of a setext heading
                const next = items[heading.item + 1];
                if (next && !next.code && /^ {0,3}(=+|-+)[ \t]*$/.test(next.line)) {
                    items[heading.item + 1] = { code: false, line: "" };
                }
            });
            const content = items
                .map((item): Item => {
                    if (item.code) {
                        return item;
                    }
                    const line = replaceLinks(item.line, (link) => {
                        if (/^[a-z][a-z0-9+.-]*:/i.test(link.target) || !link.target) {
                            return item.line.slice(link.start, link.end);
                        }
                        const bang = link.image ? "!" : "";
                        if (link.image) {
                            return `${bang}[${link.text}](${new URL(link.target, docsRawUrl + doc.file).href}${link.title})`;
                        }
                        const url = new URL(link.target, docsBlobUrl + doc.file);
                        const target = docs.findIndex((d) => docsBlobUrl + d.file === url.origin + url.pathname);
                        const anchor = decodeURIComponent(url.hash.slice(1)).toLowerCase();
                        const href = target < 0 ? url.href : placeholder(target, anchor ? `a-${anchor}` : "top");
                        return `[${link.text}](${href}${link.title})`;
                    });
                    // "section 4.2" in the text
                    return {
                        code: false,
                        line: line.replace(/\b(sections?)[ \t]+(\d+(?:\.\d+)*)\b(?![^[]*\])/gi, (match, word, old) => {
                            const target = sections.get(old);
                            return target
                                ? `${word} [${target.title}](${placeholder(index, `h-${target.heading}`)})`
                                : match;
                        }),
                    };
                })
                .map((item) => (item.code ? item.lines.join("\n") : item.line))
                .join("\n");
            return `${marker(index)}\n\n${content}`;
        })
        .join("\n\n");
}

const tags = (html: string) => html.replace(/<[^>]+>/g, "");

// Resolve the links of joinNfCoreUtilsDocs to the ids of the rendered headings, and list the docs and their sections
// for the table of contents.
export function resolveNfCoreUtilsHtml(html: string, docs: NfCoreUtilsDoc[]): { html: string; headings: TocHeading[] } {
    const parts = html.split(/<!-- nf-core-utils doc (\d+) -->/);
    const rendered = new Map<number, { depth: number; id: string; html: string }[]>();
    for (let i = 1; i < parts.length; i += 2) {
        rendered.set(
            Number(parts[i]),
            [...parts[i + 1].matchAll(/<h([1-6])\b[^>]*?\bid="([^"]*)"[^>]*>([\s\S]*?)<\/h\1>/g)].map((match) => ({
                depth: Number(match[1]),
                id: match[2],
                html: match[3],
            })),
        );
    }
    // per doc, the ids of the rendered headings by GitHub anchor of the heading in the original doc
    const anchors = docs.map((doc, index) => {
        const ids = rendered.get(index) ?? [];
        const original = headings(splitCode(doc.content));
        const slugger = new GithubSlugger();
        const byAnchor = new Map<string, string>();
        if (original.length === ids.length) {
            original.forEach((heading, number) => byAnchor.set(slugger.slug(plainText(heading.text)), ids[number].id));
        } else {
            console.warn(
                `nf-core-utils docs: ${doc.file} has ${original.length} headings but ${ids.length} once rendered, links to its sections may point to GitHub`,
            );
        }
        return { ids, byAnchor };
    });
    const resolved = html.replace(/href="#nf-core-utils-(\d+)-(top|a-[^"]*|h-\d+)"/g, (match, index, anchor) => {
        const doc = anchors[Number(index)];
        const id =
            anchor === "top"
                ? doc.ids[0]?.id
                : anchor.startsWith("h-")
                  ? doc.ids[Number(anchor.slice(2))]?.id
                  : (doc.byAnchor.get(anchor.slice(2)) ??
                    doc.ids.find((heading) => heading.id === anchor.slice(2))?.id);
        if (id) {
            return `href="#${id}"`;
        }
        console.warn(`nf-core-utils docs: no heading for the link to ${docs[Number(index)].file}#${anchor.slice(2)}`);
        return `href="${docsBlobUrl}${docs[Number(index)].file}${anchor.startsWith("a-") ? `#${anchor.slice(2)}` : ""}"`;
    });
    const toc = [...rendered.keys()]
        .sort((a, b) => a - b)
        .flatMap((index) => rendered.get(index)!)
        .filter((heading) => heading.depth <= 2)
        .map((heading) => {
            // without the link icon added by the heading anchors plugin
            const inner = heading.html.replace(/<a href="#[^"]*">\s*<i [^>]*><\/i>\s*<\/a>\s*$/, "").trim();
            const code = inner.match(/^<code\b[^>]*>([\s\S]*)<\/code>$/);
            return {
                // code only for headings that are code, as text and code don't wrap together
                text: code && !code[1].includes("</code>") ? `<code>${tags(code[1])}</code>` : tags(inner),
                slug: heading.id,
                depth: heading.depth,
            };
        });
    return { html: resolved.replace(/<!-- nf-core-utils doc \d+ -->\n?/g, ""), headings: toc };
}
