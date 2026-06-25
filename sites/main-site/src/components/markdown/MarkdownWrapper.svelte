<script lang="ts">
    import { currentHeading, Checkboxes } from "@components/store";
    import * as icons from "file-icons-js";
    import "file-icons-js/css/style.css";
    import { onMount } from "svelte";

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

        // Add file icon to Expressive Code editor frame titles
        document
            .querySelectorAll(".expressive-code figure.has-title:not(.is-terminal) figcaption .title")
            .forEach((block) => {
                // Skip blocks that already got a build-time icon (e.g. the
                // Nextflow logo injected by the file-icons EC plugin).
                if (block.querySelector("svg")) return;
                const title = block.textContent;
                const fileIcon = icons.getClass(title);
                let icon: HTMLElement;
                if (fileIcon) {
                    icon = document.createElement("span");
                    icon.classList.add("ms-1", "me-2", fileIcon);
                } else {
                    icon = document.createElement("i");
                    icon.classList.add("fa", "fa-regular", "fa-file-code", "ms-1", "me-2");
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
