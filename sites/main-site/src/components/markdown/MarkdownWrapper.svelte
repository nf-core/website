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
