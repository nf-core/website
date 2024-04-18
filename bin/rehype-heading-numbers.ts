import { visit } from 'unist-util-visit';
import type { Element } from 'hast';
import type { RehypePlugin } from '@astrojs/markdown-remark';

export const rehypeHeadingNumbers: RehypePlugin<[]> = () => {
    return (tree: any, file: any) => {
        let counters: number[] = [1];
        let minDepth = 1;

        if (file.data.astro.frontmatter.addNumbersToHeadings) {
            // set mindepth to lowest heading depth
            visit(tree, 'element', (node: Element) => {
                if (
                    node.tagName === 'h1' ||
                    node.tagName === 'h2' ||
                    node.tagName === 'h3' ||
                    node.tagName === 'h4' ||
                    node.tagName === 'h5' ||
                    node.tagName === 'h6'
                ) {
                    minDepth = Math.min(minDepth, parseInt(node.tagName.slice(1)));
                }
            });
            visit(tree, 'element', (node: Element) => {
                if (
                    node.tagName === 'h1' ||
                    node.tagName === 'h2' ||
                    node.tagName === 'h3' ||
                    node.tagName === 'h4' ||
                    node.tagName === 'h5' ||
                    node.tagName === 'h6'
                ) {
                    const depth = parseInt(node.tagName.slice(1)) - minDepth + 1;

                    // Initialize counters for new depths
                    if (depth > counters.length) {
                        // Initialize counter for new depth
                        counters.push(1);
                    } else {
                        // Increment counter for current depth and reset counters for bigger depths
                        counters = counters.slice(0, depth);
                        counters[depth - 1] = counters[depth - 1] + 1;
                    }

                    // Add counter to heading
                    const counterText = counters.slice(1, depth).join('.');

                    // remove hard-coded heading numbers
                    if (node.children[0].type === 'text') {
                        node.children[0].value = node.children[0].value.replace(/^\s*\d+(\.\d+)*\s/, '');
                    }

                    node.children.unshift({
                        type: 'text',
                        value: `${counterText} `,
                    });
                }
            });
        }
    };
};

export default rehypeHeadingNumbers;
