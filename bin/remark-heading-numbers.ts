import { visit } from 'unist-util-visit';
import type { Root, Heading } from 'mdast';
import type { RemarkPlugin } from '@astrojs/markdown-remark';

export const remarkHeadingNumbers: RemarkPlugin<[]> = () => {
    return (tree: Root, file: any) => {
        let counters: number[] = [1];
        let minDepth = 1;

        if (file.data.astro.frontmatter.addNumbersToHeadings) {
            // set mindepth to lowest heading depth
            visit(tree, 'heading', (node: Heading) => {
                minDepth = Math.min(minDepth, node.depth);
            });
            visit(tree, 'heading', (node: Heading) => {
                const depth = node.depth - minDepth + 1;

                // Initialize counters for new depths
                if (depth > counters.length) {
                    // Initialize counter for new depth
                    counters.push(1);
                } else {
                    // console.log('counters', counters, 'depth', depth, 'counters[depth-1]', counters[depth - 1]);
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
            });
        }
    };
};

export default remarkHeadingNumbers;
