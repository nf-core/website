import { visit } from 'unist-util-visit';
import type { Root } from 'mdast';
import type { RemarkPlugin } from '@astrojs/markdown-remark';

export const remarkcheckboxParser: RemarkPlugin<[]> = () => (tree: Root, file: any) => {
    let id = 0;
    let lastVisitedHeading = null;
    if (file.data.astro.frontmatter.type === 'checklist') {
        file.data.astro.frontmatter.checkboxes = [];
        visit(tree, (node: any, index: number, parent: any) => {
            if (node.type === 'heading') {
                lastVisitedHeading = node;
            }
            if (node.type === 'listItem' && node.checked !== null) {
                const uniqueId = `checkbox-${id++}`;
                // convert it to an html input element with label, with a special handling for nested lists

                if (node.children[1]?.type === 'list') {
                    node.children.unshift({
                        type: 'html',
                        value: `<div class="form-check"><label class="form-check-label"><input type="checkbox" class="form-check-input" id="${uniqueId}" name="${uniqueId}" ${node.checked ? 'checked' : ''}/> `,
                    });

                    // insert closing label after the first element and then add the nodeList
                    node.children.splice(2, 0, {
                        type: 'html',
                        value: `</label></div>`,
                    });
                } else {
                    node.children.unshift({
                        type: 'html',
                        value: `<div class="form-check"><label class="form-check-label"><input type="checkbox" class="form-check-input" id="${uniqueId}" name="${uniqueId}" ${
                            node.checked ? 'checked' : ''
                        }/> `,
                    });

                    node.children.push({
                        type: 'html',
                        value: `</label></div>`,
                    });
                }
                node.checked = null;
                file.data.astro.frontmatter.checkboxes.push({
                    id: uniqueId,
                    checked: node.checked,
                    label: node.children[1].children,
                    heading: lastVisitedHeading,
                });
            }
        });
    }
};
