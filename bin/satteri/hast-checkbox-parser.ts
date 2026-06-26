// Satteri replacement for bin/rehype-checkbox-parser.ts.
// Turns GFM task-list checkboxes into Bootstrap form-checks and collects them for
// the page (the old plugin wrote them to file.data.astro.frontmatter.checkboxes;
// here the caller reads them via getCheckboxes() after rendering).
// Only included in the pipeline of documents that opt in via markdownPlugin:
// checklist — see markdownConfig.ts.
import { h } from "hastscript";
import { defineHastPlugin } from "satteri";

export interface CheckboxEntry {
    id: string;
    checked: boolean | undefined;
    headingId: string | undefined;
}

const isElement = (node: any, tagName?: string) =>
    node?.type === "element" && (tagName === undefined || node.tagName === tagName);

export function createCheckboxParser() {
    let checkboxes: CheckboxEntry[] = [];

    const plugin = () => {
        checkboxes = [];
        let id = 0;
        let lastHeadingId: string | undefined;

        const makeInput = (input: any, uniqueId: string) => ({
            ...input,
            properties: {
                ...input.properties,
                id: uniqueId,
                name: uniqueId,
                type: "checkbox",
                className: "form-check-input",
                disabled: false,
            },
        });

        const formCheck = (children: any[]) =>
            h("div", { class: "form-check" }, [h("label", { class: "form-check-label" }, children)]);

        return defineHastPlugin({
            name: "checkbox-parser",
            element: {
                filter: ["h1", "h2", "h3", "h4", "h5", "h6", "li", "p"],
                visit(node: any) {
                    if (/^h[1-6]$/.test(node.tagName)) {
                        lastHeadingId = node.properties?.id;
                        return;
                    }

                    if (!node.children) return;

                    if (node.tagName === "li" && isElement(node.children[0], "input")) {
                        const uniqueId = `checkbox-${lastHeadingId}-${id++}`;
                        const input = makeInput(node.children[0], uniqueId);
                        const oldChildren = [input, ...node.children.slice(1)];

                        const ulIndex = node.children.findIndex((c: any) => isElement(c, "ul"));
                        const newChildren =
                            ulIndex === -1
                                ? [formCheck(oldChildren)]
                                : [formCheck(oldChildren.slice(0, ulIndex)), node.children[ulIndex]];

                        checkboxes.push({
                            id: uniqueId,
                            checked: input.properties.checked,
                            headingId: lastHeadingId,
                        });

                        return { ...node, children: newChildren };
                    }

                    // Handle checkboxes in paragraphs - look for input anywhere in children
                    if (node.tagName === "p") {
                        const inputIndex = node.children.findIndex(
                            (c: any) => isElement(c, "input") && c.properties?.type === "checkbox",
                        );
                        if (inputIndex === -1) return;

                        const uniqueId = `checkbox-${lastHeadingId}-${id++}`;
                        const input = makeInput(node.children[inputIndex], uniqueId);
                        const oldChildren = [...node.children];
                        oldChildren[inputIndex] = input;

                        checkboxes.push({
                            id: uniqueId,
                            checked: input.properties.checked,
                            headingId: lastHeadingId,
                        });

                        // Transform the paragraph into a proper Bootstrap form-check div
                        return {
                            ...node,
                            tagName: "div",
                            properties: { className: ["form-check"] },
                            children: [h("label", { class: "form-check-label" }, oldChildren)],
                        };
                    }
                },
            },
        });
    };

    return { plugin, getCheckboxes: () => checkboxes };
}
