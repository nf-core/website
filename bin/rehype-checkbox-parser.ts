import { visit } from 'unist-util-visit';
import type { RehypePlugin } from '@astrojs/markdown-remark';

export const rehypeCheckboxParser: RehypePlugin<[]> = () => (tree: any, file: any) => {
    let id = 0;
    let lastVisitedHeading = null;
    if (file.data.astro.frontmatter.markdownPlugin?.includes('checklist')) {
        file.data.astro.frontmatter.checkboxes = [];
        visit(tree, (node: any) => {
            if (node.type === 'element' && node.tagName.startsWith('h')) {
                lastVisitedHeading = node;
            }
            if (node.type === 'element' && node.tagName === 'li' && node.children[0].tagName === 'input') {
                const uniqueId = `checkbox-${lastVisitedHeading?.properties.id}-${id++}`;

                node.children[0].properties.id = uniqueId;
                node.children[0].properties.name = uniqueId;
                node.children[0].properties.type = 'checkbox';
                node.children[0].properties.class = 'form-check-input';
                node.children[0].properties.disabled = false;

                const oldChildren = node.children;

                const ulIndex = node.children.findIndex((c) => c.tagName === 'ul');

                if (ulIndex === -1) {
                    const div = {
                        type: 'element',
                        tagName: 'div',
                        properties: { class: 'form-check' },
                        children: [
                            {
                                type: 'element',
                                tagName: 'label',
                                properties: { class: 'form-check-label' },
                                children: oldChildren,
                            },
                        ],
                    };

                    node.children = [div];
                } else {
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
                                children: divChildren,
                            },
                        ],
                    };
                    node.children = [];
                    node.children.push(div);
                    node.children.push(ul);
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
