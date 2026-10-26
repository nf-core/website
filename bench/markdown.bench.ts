// Benchmarks for the shared satteri markdown pipeline in bin/satteri, which renders
// every docs page, blog post and pipeline README on the website.
// Run locally with `npm run bench`; CI runs them through CodSpeed.
import { readFileSync } from "node:fs";
import { createSatteriMarkdownProcessor } from "@astrojs/markdown-satteri";
import { bench, describe } from "vitest";
import { createSatteriPluginSets, renderMarkdown } from "../bin/satteri/markdownConfig.ts";
import { createGitHubMarkdownPlugin, trimGitHubReadme } from "../bin/satteri/mdast-github-markdown.ts";

const fixture = (name: string) => readFileSync(new URL(`./fixtures/${name}`, import.meta.url), "utf8");

const moduleGuidelines = fixture("docs-module-guidelines.md");
const writingModulesTutorial = fixture("docs-writing-modules-tutorial.md");
const pipelineReadme = fixture("pipeline-readme.md");

const shortDescription =
    "Assembly, binning and annotation of metagenomes :microbe: using [MEGAHIT](https://github.com/voutcn/megahit), `SPAdes` and **MetaBAT2**.";

const admonitionHeavy = Array.from(
    { length: 20 },
    (_, i) => `## Section ${i}

:::note{title="Note ${i}"}
Some _emphasised_ text with \`inline code{:groovy}\` and a [link](../usage.md#section-${i}).
:::

:::warning
Watch out for **this** :warning:.
:::

- [ ] Unchecked item ${i}
- [x] Checked item ${i}
`,
).join("\n");

// Mirrors renderPipelineMarkdown in sites/main-site/utils/pipelineMarkdown.ts.
async function renderPipelineMarkdown(body: string, parentDirectory: string) {
    const { features, mdastPlugins, hastPlugins } = createSatteriPluginSets();
    const processor = await createSatteriMarkdownProcessor({
        syntaxHighlight: false,
        features,
        mdastPlugins: [
            createGitHubMarkdownPlugin({ repo: "example", ref: "1.0.0", version: "1.0.0", parentDirectory }),
            ...mdastPlugins,
        ],
        hastPlugins,
    } as any);
    return processor.render(trimGitHubReadme(body));
}

describe("renderMarkdown", () => {
    bench("short pipeline description", async () => {
        await renderMarkdown(shortDescription);
    });

    bench("admonitions and checklists", async () => {
        await renderMarkdown(admonitionHeavy);
    });

    bench("docs page: module guidelines", async () => {
        await renderMarkdown(moduleGuidelines);
    });

    bench("docs page: writing modules tutorial", async () => {
        await renderMarkdown(writingModulesTutorial);
    });
});

describe("renderPipelineMarkdown", () => {
    bench("pipeline README", async () => {
        await renderPipelineMarkdown(pipelineReadme, "");
    });
});

describe("mdast-github-markdown", () => {
    bench("trimGitHubReadme", () => {
        trimGitHubReadme(pipelineReadme);
    });

    const ctx = { setProperty: (node: any, key: string, value: any) => (node[key] = value) };
    const links = [
        "usage.md#a-note",
        "../usage.md",
        "subdir/page.mdx",
        "../CHANGELOG.md#v250",
        "CITATIONS.md#bowtie2",
        "../assets/foo.csv",
        "https://nf-co.re/usage.md",
        "#anchor-only",
    ];

    bench("rewrite links", () => {
        const plugin: any = createGitHubMarkdownPlugin({
            repo: "mag",
            ref: "2.5.4",
            version: "dev",
            parentDirectory: "docs",
        });
        for (let i = 0; i < 25; i++) {
            for (const url of links) {
                plugin.link({ url }, ctx);
            }
        }
    });
});
