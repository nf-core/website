<script lang="ts">
    import { currentHeading } from "@components/store";
    import SchemaListingTocButtons from "@components/schema/SchemaListingTocButtons.svelte";
    import { onMount } from "svelte";
    interface Props {
        headings: {
            text: string;
            slug: string;
            depth: number;
            fa_icon?: string;
        }[];
        showHiddenBtn: boolean;
    }

    let { headings, showHiddenBtn }: Props = $props();
    // indent the headings by level, as in the sidebar
    const minDepth = $derived(Math.min(...headings.map((heading) => heading.depth)));

    // Bootstrap doesn't position the dropdowns of a navbar: shift the menu left when it would go past the screen
    let dropdown: HTMLElement;
    onMount(() => {
        const keepOnScreen = () => {
            const menu = dropdown.querySelector<HTMLElement>(".dropdown-menu");
            if (!menu) return;
            menu.style.left = "";
            const overflow = menu.getBoundingClientRect().right - (window.innerWidth - 8);
            if (overflow > 0) {
                menu.style.left = `${-overflow}px`;
            }
        };
        dropdown.addEventListener("shown.bs.dropdown", keepOnScreen);
        return () => dropdown.removeEventListener("shown.bs.dropdown", keepOnScreen);
    });
</script>

<!-- below xl, where the table of contents of the sidebar is hidden -->
<div class="d-xl-none toc-md">
    <div class="dropdown" bind:this={dropdown}>
        <button
            class="btn btn-sm btn-outline-secondary dropdown-toggle text-body-secondary"
            type="button"
            id="dropdownMenuButton"
            data-bs-toggle="dropdown"
            aria-expanded="false"
        >
            <i class="fa-solid fa-list" aria-hidden="true"></i> On this page
        </button>
        <ul class="dropdown-menu" aria-labelledby="dropdownMenuButton">
            {#each headings as heading (heading)}
                <li class:active={heading.slug === $currentHeading}>
                    <a
                        class="dropdown-item"
                        href={"#" + heading.slug}
                        style:padding-left={`${heading.depth - minDepth + 1}rem`}
                    >
                        {#if heading.fa_icon}
                            <i class={heading.fa_icon} aria-hidden="true"></i>
                        {/if}
                        {@html heading.text}
                    </a>
                </li>
            {/each}
            {#if showHiddenBtn}
                <SchemaListingTocButtons />
            {/if}
        </ul>
    </div>
</div>

<style lang="scss">
    .toc-md {
        .dropdown-menu {
            max-height: 90dvh;
            // as wide as the longest heading, but not wider than the screen
            width: max-content;
            max-width: calc(100vw - 1rem);
            overflow-y: auto;
        }
        // long headings wrap instead of going past the screen
        .dropdown-item {
            white-space: normal;
        }
    }
    :global(.navbar:has(.docs-nav)) .toc-md .dropdown-menu {
        z-index: 1050;
    }
    li.active {
        background-color: rgba(var(--bs-success), 0.75);
    }

    :global([data-bs-theme="dark"]) {
        li.active {
            background-color: rgba(var(--bs-success), 0.75);
        }
    }
</style>
