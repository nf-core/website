<script>
    import Markdown from 'svelte-exmarkdown';
    import emoji from 'remark-emoji';
    import remarkGfm from 'remark-gfm';
    import remarkDirective from 'remark-directive';
    import admonitionPlugin from '/bin/remark-admonition.js';
    import addClasses from 'rehype-add-classes';
    import rehypeAutolinkHeadings from 'rehype-autolink-headings';
    // import rehypePrettyCode from 'rehype-pretty-code';
    import rehypeSlug from 'rehype-slug';
    import urls from 'rehype-urls';
    import rehypeWrap from 'rehype-wrap-all';
    import { h } from 'hastscript';

    export let md;

    // replace > ⚠️ with :::warning amd add ::: at the end
    // md = md.replace(/(> ⚠️)(.*)(\n|\.$)/g, ':::note \n $2 \n:::\n');
    if (md) {
        // replace newline with <br>
        md = md.replace(/(\n)/g, '  \n');
        // escape *
        md = md.replace(/(\*.)/g, '\\$1');
    }
</script>

<Markdown
    {md}
    plugins={[
        {
            remarkPlugin: [emoji, remarkGfm, remarkDirective, admonitionPlugin],
            rehypePlugin: [
                rehypeSlug,
                [
                    rehypeAutolinkHeadings,
                    {
                        behavior: 'append',
                        content: h('i.ms-1.fas.fa-link.invisible'),
                    },
                ],
                [
                    addClasses,
                    {
                        table: 'table table-hover table-sm small',
                    },
                ],
                [
                    rehypeWrap,
                    {
                        selector: 'table',
                        wrapper: 'div.table-responsive',
                    },
                ],
                [
                    urls,
                    (url) => {
                        if (url.href?.endsWith('.md')) {
                            url.href = url.href.replace(/\.md$/, '/');
                            url.pathname = url.pathname.replace(/\.md$/, '/');
                            url.path = url.path.replace(/\.md$/, '/');
                        }
                    },
                ],
                // [ // vite doesn't like to compile rehype-pretty-code
                //     rehypePrettyCode,
                //     {
                //         langPrefix: 'language-',
                //     },
                // ],
            ],
        },
    ]}
/>
