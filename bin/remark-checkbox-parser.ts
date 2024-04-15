import { visit } from 'unist-util-visit';
import type { Root } from 'mdast';
import type { RemarkPlugin } from '@astrojs/markdown-remark';

let id = 0;

export const remarkcheckboxParser: RemarkPlugin<[]> = () => (tree: Root, file: any) => {
    if (file.data.astro.frontmatter.type === 'checklist') {
        file.data.astro.frontmatter.checkboxes = [];
        visit(tree, 'listItem', (node: any) => {
            if (node.checked !== null) {
                console.log(node);
                const uniqueId = `checkbox-${id++}`;
                // convert it to an html input element with label
                node.children.unshift({
                    type: 'html',
                    value: `<input type="checkbox" id="${uniqueId}" name="${uniqueId}" ${
                        node.checked ? 'checked' : ''
                    }/> <label for="${uniqueId}">`,
                });

                node.children.push({
                    type: 'html',
                    value: `</label>`,
                });
                node.checked = null;
                file.data.astro.frontmatter.checkboxes.push({
                    id: uniqueId,
                    checked: node.checked,
                    label: node.children[1].children,
                });
            }
        });
    }
};
