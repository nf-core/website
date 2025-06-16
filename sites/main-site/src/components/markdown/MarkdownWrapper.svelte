<script lang="ts">
    import { currentHeading, Checkboxes } from "@components/store";
    import * as icons from "file-icons-js";
    import "file-icons-js/css/style.css";
    import { onMount, mount } from "svelte";
    import CopyButton from "@components/CopyButton.svelte";

    interface Props {
        headings?: { text: string; slug: string; depth: number; fa_icon?: string }[];
        children?: import("svelte").Snippet;
    }
    let { headings = [], children }: Props = $props();
    onMount(() => {
        if (typeof window !== "undefined" && window.location.pathname.includes("/events/") && !window.location.hash) {
            window.scrollTo(0, 0);
        }
        async function renderDiagrams(graphs, mermaidModule) {
            mermaidModule.initialize({
                startOnLoad: false,
                fontFamily: "var(--sans-font)",
                // @ts-ignore This works, but TS expects a enum for some reason
                theme: document.documentElement.getAttribute("data-bs-theme") === "dark" ? "dark" : "neutral",
            });

            for (const graph of graphs) {
                const content = graph.getAttribute("data-content");
                if (!content) continue;
                let svg = document.createElement("svg");
                const id = (svg.id = "mermaid-" + Math.round(Math.random() * 100000));
                graph.appendChild(svg);
                mermaidModule.render(id, content).then((result) => {
                    graph.innerHTML = result.svg;
                });
            }
        }
        const graphs = document.getElementsByClassName("mermaid");
        if (document.getElementsByClassName("mermaid").length > 0) {
            // Dynamically import mermaid only when needed
            import("mermaid").then((mermaidModule) => {
                renderDiagrams(graphs, mermaidModule.default);
                window.addEventListener("theme-changed", () => {
                    renderDiagrams(graphs, mermaidModule.default);
                });
            });
        }

        // find current heading in viewport with IntersectionObserver
        const observer = new IntersectionObserver(
            (entries) => {
                entries.forEach((entry) => {
                    if (entry.isIntersecting) {
                        currentHeading.set(entry.target.id);
                    }
                });
            },
            {
                rootMargin: "0px 0px -92% 0px",
            },
        );
        headings.forEach((heading) => {
            const element = document.querySelector("#" + heading.slug);
            if (element) {
                observer.observe(element);
            }
        });

        // Add "Copy code" button in code blocks
        const copyButtonLabel = "<i class='fa-regular fa-clipboard'></i>";
        const copiedButtonLabel = `<span class='font-sans-serif'>Copied </span><i class='fa-regular fa-clipboard-check'></i>`;
        document
            .querySelectorAll(
                "figure[data-rehype-pretty-code-figure] pre:not([data-language='console']):not([data-language='tree']):not([data-language='plaintext'])",
            )
            .forEach((block) => {
                block.classList.add("position-relative");
                const copyText = block.querySelector("code")?.innerText;
                if (copyText) {
                    // check if block has only one child, i.e. is a single line code block, so we need less top and bottom margin for button
                    const SingleLine = block.children[0].childElementCount === 1 ? "single-line" : "";
                    mount(CopyButton, {
                        target: block, // Specify the target element for the Svelte component
                        props: {
                            text: copyText,
                            label: copyButtonLabel,
                            copiedLabel: copiedButtonLabel,
                            classes:
                                SingleLine +
                                " copy-code-button btn btn-sm btn-outline-secondary position-absolute top-0 end-0 opacity-50",
                            copiedClasses:
                                SingleLine +
                                " copy-code-button btn btn-sm btn-outline-success position-absolute top-0 end-0",
                        },
                    });
                }
            });
        // Add file icon to code block titles
        document.querySelectorAll("[data-rehype-pretty-code-title]").forEach((block) => {
            const title = block.textContent;
            const fileIcon = icons.getClass(title);
            let icon: HTMLElement;
            if (fileIcon) {
                icon = document.createElement("span");
                icon.classList.add("ms-1", "me-2", fileIcon);
            } else {
                icon = document.createElement("i");
                icon.classList.add("fa-regular", "fa-file-code", "ms-1", "me-2");
            }
            block.prepend(icon);
        });

        // change stored Checkboxes to checked
        $Checkboxes.forEach((checkbox) => {
            const element = document.getElementById(checkbox.id);
            if (element) {
                (element as HTMLInputElement).checked = true;
            }
        });

        // Update Checkboxes store when checkboxes are clicked
        const checkboxes = document.querySelectorAll('input[type="checkbox"]');
        checkboxes.forEach((checkbox) => {
            // First, set initial state from store, but only if not all checkboxes would be checked
            if (!$Checkboxes.every((c) => c.checked)) {
                const storedCheckbox = $Checkboxes.find((c) => c.id === checkbox.id);
                if (storedCheckbox) {
                    (checkbox as HTMLInputElement).checked = true;
                }
            }
            checkbox.addEventListener("change", () => {
                const checkedBoxes = Array.from(checkboxes)
                    .filter((cb) => (cb as HTMLInputElement).checked)
                    .map((cb) => ({
                        id: cb.id,
                        checked: true,
                    }));

                Checkboxes.set(checkedBoxes);
            });
        });
    });
</script>

<div>
    {@render children?.()}
</div>

<style lang="scss">
    :global(.file-tree) {
        // CSS Custom Properties with modern light-dark() function
        --tree-x-space: 1.5rem;
        --tree-y-space: 0.125rem;
        --tree-y-pad: 0;
        --tree-border-color: light-dark(var(--bs-gray-300), var(--bs-gray-600));
        --tree-text-muted: light-dark(var(--bs-gray-600), var(--bs-gray-400));
        --tree-text-subtle: light-dark(var(--bs-gray-500), var(--bs-gray-400));

        font-size: 0.875rem;
        overflow-x: auto;
        container-type: inline-size;

        // Nested list styling with modern CSS nesting
        :global(.nested-list) {
            border-inline-start: 1px solid var(--tree-border-color);
            padding-block: 0;
            padding-inline: 0 0 0 0.5rem;
        }

        // Tree entry components
        :global(.tree-entry) {
            display: flex;
            align-items: center;
            gap: 0.25rem;
        }

        :global(.file-name) {
            font-weight: 500;
            color: var(--bs-body-color);
        }

        :global(.comment) {
            color: var(--tree-text-muted);
            font-style: italic;
            margin-inline-start: auto;

            // Container query for responsive layout
            @container (max-width: 30rem) {
                display: block;
                margin-inline-start: 1.5rem;
                margin-block-start: 0.25rem;
            }
        }

        // Icon styling with logical properties
        :global(.tree-icon) {
            inline-size: 0.875rem;
            block-size: 0.875rem;
            flex-shrink: 0;

            :global(path) {
                fill: var(--tree-text-muted) !important;
            }
        }

        // Highlighted entries
        :global(.highlight) {
            display: inline-block;
            border-radius: 0.25rem;
            padding-inline-end: 0.5rem;
            color: white;
            background-color: var(--bs-primary);

            .tree-icon path {
                fill: currentColor !important;
            }
        }

        // Empty state
        :global(.empty) {
            color: var(--tree-text-subtle);
            padding-inline-start: 0.5rem;
            font-style: italic;
        }

        // Tree item spacing
        :global(.tree-item) {
            margin-block: var(--tree-y-space);
            padding-left: 0.5rem;

            :global(> div) {
                padding-block: 0.15rem 0.5rem;
            }
        }
    }

    // Fallback for browsers without light-dark() support
    @supports not (color: light-dark(white, black)) {
        :global(.file-tree) {
            --tree-border-color: var(--bs-gray-300);
            --tree-text-muted: var(--bs-gray-600);
            --tree-text-subtle: var(--bs-gray-500);
        }

        :global([data-bs-theme="dark"] .file-tree) {
            --tree-border-color: var(--bs-gray-600);
            --tree-text-muted: var(--bs-gray-400);
            --tree-text-subtle: var(--bs-gray-400);
        }
    }
</style>
