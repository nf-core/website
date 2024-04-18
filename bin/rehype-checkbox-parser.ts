import { visit } from 'unist-util-visit';
import type { RehypePlugin } from '@astrojs/markdown-remark';

export const rehypeCheckboxParser: RehypePlugin<[]> = () => (tree: any, file: any) => {
    let id = 0;
    let lastVisitedHeading = null;
    if (file.data.astro.frontmatter.type === 'checklist') {
        file.data.astro.frontmatter.checkboxes = [];
        visit(tree, (node: any) => {
            if (node.type === 'element' && node.tagName.startsWith('h')) {
                lastVisitedHeading = node;
            }
            if (node.type === 'element' && node.tagName === 'li') {
                const uniqueId = `checkbox-${lastVisitedHeading?.properties.id}-${id++}`;
                const oldChildren = node.children.filter((c) => c.tagName !== 'input');
                // convert it to an bootstrap checkbox element with label, with a special handling for nested lists

                const ulIndex = node.children.findIndex((c) => c.tagName === 'ul');
                if (ulIndex !== -1) {
                    // wrap all children except the ul in a div
                    const ul = node.children[ulIndex];
                    const divChildren = oldChildren.slice(0, ulIndex);
                    const div = {
                        type: 'element',
                        tagName: 'div',
                        properties: { class: 'form-check' },
                        children: [
                            {
                                type: 'element',
                                tagName: 'label',
                                properties: { class: 'form-check-label' },
                                children: [
                                    {
                                        type: 'element',
                                        tagName: 'input',
                                        properties: {
                                            type: 'checkbox',
                                            class: 'form-check-input',
                                            id: uniqueId,
                                            name: uniqueId,
                                        },
                                    },
                                ],
                            },
                            ...divChildren,
                        ],
                    };
                    node.children = [];
                    node.children.push(div);
                    node.children.push(ul);
                } else {
                    node.children = [
                        {
                            type: 'element',
                            tagName: 'div',
                            properties: { class: 'form-check' },
                            children: [
                                {
                                    type: 'element',
                                    tagName: 'label',
                                    properties: { class: 'form-check-label' },
                                    children: [
                                        {
                                            type: 'element',
                                            tagName: 'input',
                                            properties: {
                                                type: 'checkbox',
                                                class: 'form-check-input',
                                                id: uniqueId,
                                                name: uniqueId,
                                            },
                                        },
                                        ...oldChildren,
                                    ],
                                },
                            ],
                        },
                    ];
                }
                file.data.astro.frontmatter.checkboxes.push({
                    id: uniqueId,
                    checked: node.properties.checked,
                    headingId: lastVisitedHeading?.properties.id,
                });
            }
        });
    }
};
