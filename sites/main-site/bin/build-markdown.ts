import admonitionsPlugin from './remark-admonitions.js';
import { mermaid } from './remark-mermaid.ts';
import { unified } from 'unified';
import remarkParse from 'remark-parse';
import addClasses from 'rehype-add-classes';
import rehypeAutolinkHeadings from 'rehype-autolink-headings';
import rehypeKatex from 'rehype-katex';
import rehypePrettyCode from 'rehype-pretty-code';
import rehypeSlug from 'rehype-slug';
import urls from 'rehype-urls';
import rehypeWrap from 'rehype-wrap-all';
import remarkDirective from 'remark-directive';
import emoji from 'remark-emoji';
import remarkGfm from 'remark-gfm';
import remarkMath from 'remark-math';

// take markdown string as input and return html string as output

export async function renderMarkdown(markdown) {
    const processor = unified()
        .use(remarkParse)
        .use(remarkDirective)
        .use(remarkGfm)
        .use(remarkMath)
        .use(emoji)
        .use(admonitionsPlugin)
        .use(mermaid)
        .use(remark2rehype)
        .use(rehypeAutolinkHeadings, {
            behavior: 'wrap',
            properties: {
                className: ['anchor'],
            },
        })
        .use(rehypeSlug)
        .use(rehypeKatex)
        .use(rehypePrettyCode)
        .use(
            urls,
            (url) => {
                if (url.startsWith('/')) {
                    return url;
                }
                return `https://deno.com/${url}`;
            },
            (node) => node.tagName === 'a',
        )
        .use(addClasses, {
            pre: 'prettyprint',
        })
        .use(rehypeWrap, {
            selector: 'table',
            wrapper: (el) => {
                const wrapper = document.createElement('div');
                wrapper.classList.add('table-wrapper');
                el.before(wrapper);
                wrapper.append(el);
            },
        });

    const result = await processor.process(markdown);
    return String(result);
}
