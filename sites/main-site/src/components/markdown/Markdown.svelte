<script lang="ts">
    import Markdown from "svelte-exmarkdown";
    import emoji from "remark-emoji";
    import remarkGfm from "remark-gfm";
    import remarkDirective from "remark-directive";
    import admonitionsPlugin from "@root/bin/remark-admonitions.js";
    import addClasses from "rehype-class-names";
    import rehypeAutolinkHeadings from "rehype-autolink-headings";
    import rehypeHighlight from "rehype-highlight";
    import rehypeSlug from "rehype-slug";
    import rehypeWrap from "rehype-wrap-all";
    import { h } from "hastscript";
    import remarkMath from "remark-math";
    import rehypeKatex from "rehype-katex";
    import type { Plugin } from "svelte-exmarkdown";

    import "../../../../../node_modules/highlight.js/styles/github-dark.css";
    let { md } = $props();

    let plugins: Plugin[] = [];
    plugins.push({ remarkPlugin: [emoji] });
    plugins.push({ remarkPlugin: [remarkGfm] });
    plugins.push({ remarkPlugin: [remarkDirective] });
    plugins.push({ remarkPlugin: [admonitionsPlugin] });
    plugins.push({ remarkPlugin: [remarkMath] });
    plugins.push({ rehypePlugin: [rehypeSlug] });

    plugins.push({
        rehypePlugin: [
            rehypeAutolinkHeadings,
            {
                behavior: "append",
                content: h("i.ms-1.fas.fa-link.fa-xs.invisible"),
            },
        ],
    });
    plugins.push({
        rehypePlugin: [
            addClasses,
            {
                table: "table table-hover table-sm small",
            },
        ],
    });
    plugins.push({
        rehypePlugin: [
            rehypeWrap,
            {
                selector: "table",
                wrapper: "div.table-responsive",
            },
        ],
    });
    // plugins.push({
    //     rehypePlugin: [
    //         urls,
    //         (url) => {
    //             const regex = /^https:\/\/(raw.)*github/;
    //             if (!regex.test(url.href) && url.href?.endsWith(".md")) {
    //                 url.href = url.href.replace(/\.md$/, "/");
    //                 url.pathname = url.pathname.replace(/\.md$/, "/");
    //                 url.path = url.path.replace(/\.md$/, "/");
    //             } else if (!regex.test(url.href) && url.href?.endsWith(".mdx")) {
    //                 url.href = url.href.replace(/\.mdx$/, "/");
    //                 url.pathname = url.pathname.replace(/\.mdx$/, "/");
    //                 url.path = url.path.replace(/\.mdx$/, "/");
    //             }
    //         },
    //     ],
    // });
    plugins.push({
        rehypePlugin: [rehypeKatex],
    });
    plugins.push({
        rehypePlugin: [rehypeHighlight],
    });
</script>

<Markdown {md} {plugins} />
