// Link rewriting checks for createGitHubMarkdownPlugin: run with `npm run test-unit`.
// Pipeline docs pages are served at trailing-slash URLs, so a relative href to a
// sibling doc resolves one level too deep and 404s — see nf-core/website#4384.
import assert from "node:assert/strict";
import test from "node:test";
import { createGitHubMarkdownPlugin } from "./mdast-github-markdown.ts";

const ctx = { setProperty: (node: any, key: string, value: any) => (node[key] = value) };

function rewriteLink(url: string, parentDirectory: string): string {
    const plugin: any = createGitHubMarkdownPlugin({ repo: "mag", ref: "2.5.4", version: "dev", parentDirectory });
    const node: any = { url };
    plugin.link(node, ctx);
    return node.url;
}

test("markdown links point at the same pipeline and version", () => {
    assert.equal(rewriteLink("usage.md#a-note", "docs"), "/mag/dev/docs/usage/#a-note");
    assert.equal(rewriteLink("usage.md", "docs"), "/mag/dev/docs/usage/");
    assert.equal(rewriteLink("../usage.md", "docs/subdir"), "/mag/dev/docs/usage/");
    assert.equal(rewriteLink("subdir/page.mdx", "docs"), "/mag/dev/docs/subdir/page/");
    // README (parentDirectory "") linking into the docs pages
    assert.equal(rewriteLink("docs/usage.md", ""), "/mag/dev/docs/usage/");
    assert.equal(rewriteLink("../README.md", "docs"), "/mag/dev/");
});

test("files without a page on the website still go to GitHub", () => {
    assert.equal(rewriteLink("../CHANGELOG.md", "docs"), "https://github.com/nf-core/mag/blob/2.5.4/CHANGELOG.md");
    // anchors survive: the reader lands on the release section, not the top of the changelog
    assert.equal(
        rewriteLink("../CHANGELOG.md#v250", "docs"),
        "https://github.com/nf-core/mag/blob/2.5.4/CHANGELOG.md#v250",
    );
    assert.equal(
        rewriteLink("CITATIONS.md#bowtie2", ""),
        "https://github.com/nf-core/mag/blob/2.5.4/CITATIONS.md#bowtie2",
    );
    assert.equal(rewriteLink("CITATIONS.md", ""), "https://github.com/nf-core/mag/blob/2.5.4/CITATIONS.md");
    assert.equal(rewriteLink("../assets/foo.csv", "docs"), "https://github.com/nf-core/mag/blob/2.5.4/assets/foo.csv");
});

test("absolute and anchor-only links are left alone", () => {
    assert.equal(rewriteLink("https://nf-co.re/usage.md", "docs"), "https://nf-co.re/usage.md");
    assert.equal(rewriteLink("#anchor-only", "docs"), "#anchor-only");
});
