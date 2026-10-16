// Checks for the adaptation of the nf-core-utils docs to /nf-core-utils/: run with `npm run test-unit`.
import assert from "node:assert/strict";
import test from "node:test";
import { markdownToHtml } from "satteri";
import {
    admonitionsHastPlugin,
    admonitionsMdastPlugin,
    directiveTextMergePlugin,
} from "../../../bin/satteri/admonitions.ts";
import { headingAnchorsPlugin } from "../../../bin/satteri/hast-heading-anchors.ts";
import { markdownLinksPlugin } from "../../../bin/satteri/hast-markdown-links.ts";
import { tablesPlugin } from "../../../bin/satteri/hast-tables.ts";
import { emojiPlugin } from "../../../bin/satteri/mdast-emoji.ts";
import {
    docsBlobUrl,
    docsRawUrl,
    joinNfCoreUtilsDocs,
    resolveNfCoreUtilsHtml,
    splitCode,
    toWebsiteMarkdown,
    type NfCoreUtilsDoc,
} from "./nf-core-utils.ts";

// the page, rendered with the plugins of the website that change text and headings (bin/satteri/markdownConfig.ts
// also has code highlighting, mermaid and math, which import files node can't load)
async function page(docs: NfCoreUtilsDoc[]) {
    const { html } = await markdownToHtml(joinNfCoreUtilsDocs(docs), {
        features: { directive: true, gfm: true, smartPunctuation: true },
        mdastPlugins: [emojiPlugin, directiveTextMergePlugin, admonitionsMdastPlugin] as any[],
        hastPlugins: [admonitionsHastPlugin, headingAnchorsPlugin(), tablesPlugin, markdownLinksPlugin] as any[],
    });
    return resolveNfCoreUtilsHtml(html, docs);
}

const ids = (html: string) => [...html.matchAll(/<h[1-6][^>]*\bid="([^"]*)"/g)].map((match) => match[1]);
const hrefs = (html: string) =>
    [...html.matchAll(/<a href="([^"]*)"[^>]*>([^<]*)</g)].map((match) => [match[2], match[1]]);

test("code fences close with the same fence, and at the end of a doc", () => {
    const blocks = splitCode("````md\n```bash\n# not a heading\n```\n````\n# heading\n```groovy\nopen");
    assert.deepEqual(
        blocks.map((block) => (block.code ? block.lines.length : block.line)),
        [5, "# heading", 3],
    );
    assert.equal(splitCode("```inline``` code")[0].code, false);
    assert.equal(splitCode("```\n```groovy\n```").length, 1, "an info string doesn't close a fence");
});

test("MkDocs admonitions and GitHub alerts become website admonitions", () => {
    assert.equal(
        toWebsiteMarkdown('!!! tip "Where to Start"\nNew here?\n\nNext'),
        ':::tip{title="Where to Start"}\nNew here?\n:::\n\nNext',
    );
    assert.equal(
        toWebsiteMarkdown('!!! note "T"\n\n    Indented\n\n    more\n\nOut'),
        ':::note{title="T"}\nIndented\n\nmore\n:::\n\nOut',
    );
    assert.equal(toWebsiteMarkdown('!!! warning "T" - a - b'), ':::warning{title="T"}\n- a\n- b\n:::');
    assert.equal(toWebsiteMarkdown("!!! tip 'Single'\ntext"), ':::tip{title="Single"}\ntext\n:::');
    assert.equal(toWebsiteMarkdown('!!! bug "a {b}"\nx'), ':::danger{title="a (b)"}\nx\n:::');
    assert.equal(toWebsiteMarkdown('- item\n\n  !!! tip "T"\n  text'), '- item\n\n  :::tip{title="T"}\n  text\n  :::');
    // the code block of the paragraph stays whole, blank line included
    assert.equal(
        toWebsiteMarkdown('!!! tip "T"\ntext\n```\na\n\nb\n```\n\nNext'),
        ':::tip{title="T"}\ntext\n```\na\n\nb\n```\n:::\n\nNext',
    );
    assert.equal(
        toWebsiteMarkdown('> [!TIP] "Quoted"\n>\n> - a\n> b\nlazy\n\nOut'),
        ':::tip{title="Quoted"}\n- a\nb\nlazy\n:::\n\nOut',
    );
    assert.equal(toWebsiteMarkdown('> [!WARNING] Say "hi"\n> text'), ":::warning{title=\"Say 'hi'\"}\ntext\n:::");
    assert.equal(toWebsiteMarkdown('```\n!!! tip "T"\n```'), '```\n!!! tip "T"\n```');
});

test("emojis go from the text and code comments, not from code", () => {
    assert.equal(toWebsiteMarkdown("- **🔧 Configuration**: ok ⚠️"), "- **Configuration**: ok ");
    assert.equal(toWebsiteMarkdown("👍🏽 ok, 👩‍💻 dev"), "ok, dev");
    assert.equal(toWebsiteMarkdown("A ↔ B, ℹ, ‼, © ® ™"), "A ↔ B, ℹ, ‼, © ® ™");
    assert.equal(toWebsiteMarkdown("`✅ kept` and `` a ` ✅ `` 🚀"), "`✅ kept` and `` a ` ✅ `` ");
    assert.equal(
        toWebsiteMarkdown('```groovy\n// ✅ CORRECT\nlog.info "✅ done" # ❌\nurl = "#✅"\n```'),
        '```groovy\n// CORRECT\nlog.info "✅ done" # \nurl = "#✅"\n```',
    );
    assert.equal(
        toWebsiteMarkdown('!!! tip "T"\n```\n"✅ done" // ✅ ok\n```'),
        ':::tip{title="T"}\n```\n"✅ done" // ok\n```\n:::',
    );
});

test("the docs are numbered as chapters and sections across the page", async () => {
    const { html, headings } = await page([
        {
            file: "A.md",
            content: "# A\n\n## 1. Overview\n\n### 1.1. Detail\n\n## 2024 Release plan\n\nSetext\n------\n",
        },
        { file: "utilities/B.md", content: "# B\n\n## Overview\n\nsee section 1.1\n" },
    ]);
    assert.deepEqual(
        headings.map((heading) => `${heading.depth} ${heading.text}`),
        ["1 1. A", "2 1.1. Overview", "2 1.2. 2024 Release plan", "2 1.3. Setext", "1 2. B", "2 2.1. Overview"],
    );
    assert.equal(new Set(ids(html)).size, ids(html).length);
    assert.ok(!html.includes("<!-- nf-core-utils doc"), "no markers left");
    assert.ok(headings.every((heading) => ids(html).includes(heading.slug)));
    // "section 1.1" in B is not the 1.1 of A
    assert.ok(html.includes("see section 1.1"));
});

test("table of contents and links use the ids of the rendered headings", async () => {
    const docs = [
        {
            file: "NfCoreUtilities.md",
            content: [
                "# Guide",
                "## _Emph_ heading",
                "## a -- b",
                "## Trailing ##",
                "## :rocket: launch",
                "## 🚀 Quick Start",
                "### 4.2. Advanced",
                'Links: [e](#emph-heading), [t](#Trailing), [q](#-quick-start), [s](utilities/B.md#Second "title"),',
                "[top](utilities/B.md), [gone](#nowhere), ![img](./img.png), [x](../README.md), `[c](#emph-heading)`,",
                "the patterns shown in section 4.2.",
                "```",
                "[c](#emph-heading)",
                "```",
                "```unclosed",
            ].join("\n\n"),
        },
        { file: "utilities/B.md", content: "# B\n\n## Second\n" },
    ];
    const { html, headings } = await page(docs);
    const rendered = ids(html);
    assert.ok(
        headings.every((heading) => rendered.includes(heading.slug)),
        JSON.stringify([headings, rendered]),
    );
    assert.equal(headings.at(-1)?.text, "2.1. Second", "the unclosed fence doesn't swallow the next doc");
    const links = Object.fromEntries(hrefs(html));
    for (const name of ["e", "t", "q", "s", "top"]) {
        assert.ok(links[name]?.startsWith("#") && rendered.includes(links[name].slice(1)), `${name}: ${links[name]}`);
    }
    assert.equal(links["Advanced"], "#advanced");
    assert.equal(links["gone"], `${docsBlobUrl}NfCoreUtilities.md#nowhere`);
    assert.equal(links["x"], "https://github.com/nf-core/nf-core-utils/blob/main/README.md");
    assert.ok(html.includes(`src="${docsRawUrl}img.png"`));
    assert.equal((html.match(/\[c\]\(#emph-heading\)/g) ?? []).length, 2, "links in code stay as they are");
});
