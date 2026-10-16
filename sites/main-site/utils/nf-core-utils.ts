// Fetch the user docs of the nf-core-utils plugin for /nf-core-utils/, and adapt their MkDocs syntax to the website.

export const docs_raw_url = "https://raw.githubusercontent.com/nf-core/nf-core-utils/refs/heads/main/docs/";
export const docs_blob_url = "https://github.com/nf-core/nf-core-utils/blob/main/docs/";
export const docs_entry = "NfCoreUtilities.md";

export interface NfCoreUtilsDoc {
    file: string; // path relative to docs/
    content: string;
}

// Start from the guide and follow its links to the docs of the extensions and utilities,
// the other files in docs/ are for the developers of the plugin.
// Docs at the top of docs/ come first, then those of each subfolder.
export async function fetchNfCoreUtilsDocs(): Promise<NfCoreUtilsDoc[]> {
    const docs: NfCoreUtilsDoc[] = [];
    const queue = [docs_entry];
    while (queue.length) {
        const file = queue.shift()!;
        if (docs.some((d) => d.file === file)) {
            continue;
        }
        const content = await fetch(docs_raw_url + file)
            .then((res) => (res.ok ? res.text() : Promise.reject(`Failed to fetch: ${res.status}`)))
            .catch(() =>
                file === docs_entry ? "# Error loading content\n\nUnable to load documentation from GitHub." : null,
            );
        // a broken link keeps pointing to GitHub
        if (content === null) {
            continue;
        }
        docs.push({ file, content: toWebsiteMarkdown(content) });
        for (const [, link] of stripCodeBlocks(content).matchAll(/\]\(([^)\s#]+\.md)(#[^)\s]*)?\)/g)) {
            const path = new URL(link, docs_raw_url + file).href;
            if (path.startsWith(docs_raw_url)) {
                queue.push(path.slice(docs_raw_url.length));
            }
        }
    }
    return docs.sort((a, b) => Number(a.file.includes("/")) - Number(b.file.includes("/")));
}

export function stripCodeBlocks(md: string): string {
    return md.replace(/^ *(```|~~~)[\s\S]*?^ *\1/gm, "");
}

// MkDocs and GitHub admonition types to the website ones
const admonitionTypes: Record<string, string> = {
    note: "note",
    abstract: "note",
    info: "info",
    question: "info",
    example: "info",
    important: "info",
    tip: "tip",
    success: "tip",
    warning: "warning",
    caution: "danger",
    danger: "danger",
    failure: "danger",
    bug: "danger",
};
const indented = /^( {4}|\t)/;

// Outside of code blocks:
// - turn MkDocs admonitions into website ones (`!!! tip "Title"` → `:::tip{title="Title"}`), with as content either
//   the indented lines after it, the paragraph right after it, or a `- a - b` list on the same line
// - turn GitHub alerts (`> [!TIP]`) into website admonitions
// - remove emojis
// and in code blocks, remove the emojis marking comments
export function toWebsiteMarkdown(md: string): string {
    const lines = md.split("\n");
    const out: string[] = [];
    let fence: string | undefined;
    for (let i = 0; i < lines.length; i++) {
        const line = lines[i];
        const fence_match = line.match(/^ *(`{3,}|~{3,})/);
        if (fence) {
            out.push(stripCommentEmojis(line));
            if (fence_match && fence_match[1].startsWith(fence)) {
                fence = undefined;
            }
            continue;
        }
        if (fence_match) {
            fence = fence_match[1];
            out.push(line);
            continue;
        }
        // GitHub alert, `> [!TIP] "Title"` and the quoted lines after it
        const alert = line.match(/^> *\[!(\w+)\](?: +["“]?([^"”]*?)["”]?)? *$/);
        if (alert) {
            const body: string[] = [];
            while (i + 1 < lines.length && lines[i + 1].startsWith(">")) {
                body.push(lines[++i].replace(/^> ?/, ""));
            }
            while (body.length && !body[0].trim()) {
                body.shift();
            }
            out.push(admonitionStart(alert[1], alert[2]), ...body.map(stripEmojis), ":::");
            continue;
        }
        const admonition = line.match(/^(?:!!!|\?\?\?\+?) *(\w+)(?: +["“]([^"”]*)["”])?(.*)$/);
        if (!admonition) {
            out.push(stripEmojis(line));
            continue;
        }
        const [, type, title, rest] = admonition;
        const inline = rest.trim();
        const body = inline.startsWith("- ")
            ? inline
                  .slice(2)
                  .split(/ - (?=\S)/)
                  .map((item) => `- ${item.trim()}`)
            : inline
              ? [inline]
              : [];
        let next = i + 1;
        while (next < lines.length && !lines[next].trim()) {
            next++;
        }
        if (next < lines.length && indented.test(lines[next])) {
            while (i + 1 < lines.length && (!lines[i + 1].trim() || indented.test(lines[i + 1]))) {
                body.push(lines[++i].replace(indented, ""));
            }
        } else if (!body.length) {
            while (i + 1 < lines.length && lines[i + 1].trim()) {
                body.push(lines[++i]);
            }
        }
        while (body.length && !body[body.length - 1].trim()) {
            body.pop();
        }
        out.push(admonitionStart(type, title), ...body.map(stripEmojis), ":::");
    }
    return out.join("\n");
}

function admonitionStart(type: string, title?: string): string {
    const attributes = title?.trim()
        ? `{title="${stripEmojis(title).replace(/`/g, "").trim().replace(/"/g, "'")}"}`
        : "";
    return `:::${admonitionTypes[type.toLowerCase()] ?? "note"}${attributes}`;
}

// an emoji (with its variation selector and ZWJ sequence) and the space after it, but not the copyright, registered and trade mark signs
const emoji = String.raw`(?![\u00a9\u00ae\u2122])\p{Extended_Pictographic}\ufe0f?(?:\u200d\p{Extended_Pictographic}\ufe0f?)* ?`;

// emojis outside of inline code
function stripEmojis(line: string): string {
    return line
        .split(/(`[^`]*`)/)
        .map((part, index) => (index % 2 ? part : part.replace(new RegExp(emoji, "gu"), "")))
        .join("");
}

// emojis marking code comments, as in `// ✅ CORRECT`, the ones in strings can be part of the example
function stripCommentEmojis(line: string): string {
    return line.replace(new RegExp(String.raw`((?:\/\/|#)\s*)${emoji}`, "gu"), "$1");
}
