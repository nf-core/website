// Link rewriting checks for markdownLinksPlugin: run with `npm run test-unit`.
import assert from "node:assert/strict";
import test from "node:test";
import { markdownLinksPlugin } from "./hast-markdown-links.ts";

const ctx = { setProperty: (node: any, key: string, value: any) => (node.properties[key] = value) };

function rewriteHref(href: string): string {
    const plugin: any = markdownLinksPlugin;
    const node: any = { tagName: "a", properties: { href } };
    plugin.element.visit(node, ctx);
    return node.properties.href;
}

test("links to markdown files become directory-style URLs", () => {
    assert.equal(rewriteHref("page.md"), "page/");
    assert.equal(rewriteHref("../subdir/page.mdx"), "../subdir/page/");
});

test("anchors and queries are kept", () => {
    // nf-core/website#4484: these used to keep the .md and 404
    assert.equal(rewriteHref("./vs-code.md#nf-core-extension-pack"), "./vs-code/#nf-core-extension-pack");
    assert.equal(rewriteHref("page.md?tab=usage#a-note"), "page/?tab=usage#a-note");
});

test("other links are left alone", () => {
    assert.equal(rewriteHref("#getcitation"), "#getcitation");
    assert.equal(
        rewriteHref("https://github.com/nf-core/rnaseq/blob/3.14.0/docs/usage.md#a"),
        "https://github.com/nf-core/rnaseq/blob/3.14.0/docs/usage.md#a",
    );
    assert.equal(rewriteHref("notes.md.backup"), "notes.md.backup");
    assert.equal(rewriteHref("page?file=notes.md"), "page?file=notes.md");
    assert.equal(rewriteHref("#notes.md"), "#notes.md");
    assert.equal(rewriteHref("https://example.com/readme.mdown"), "https://example.com/readme.mdown");
});
