<script lang="ts">
    import { currentHeading, showHidden, Checkboxes } from "@components/store";
    import { onMount } from "svelte";
    import ProgressIndicator from "@components/sidebar/ProgressIndicator.svelte";

    // Add snippet props to replace slots
    export let rightSidebarTop: () => any = () => null;
    export let rightSidebarLinkBar: () => any = () => null;
    export let defaultContent: () => any = () => null;

    export let headings: {
        text: string;
        slug: string;
        depth: number;
        fa_icon?: string;
        hidden?: boolean;
        checkboxes?: {
            id: string;
            label: string;
            checked: boolean;
        }[];
    }[];
    export let minNumHeadings: number = 2;
    export let minHeadingDepth: number = 1;
    export let maxHeadingDepth: number = 4;

    // filter out headings that are lower than min_heading_depth
    headings = headings.filter((h) => h.depth >= minHeadingDepth);

    // filter out headings that are higher than max_heading_depth
    headings = headings.filter((h) => h.depth <= maxHeadingDepth);

    minHeadingDepth = Math.min(...headings.map((h) => h.depth));

    const showToc = headings.length > minNumHeadings || headings.some((h) => h.checkboxes);

    let activeHeading = {};
    let hCheckboxes = headings
        .map((h) => h.checkboxes)
        .map((c) => c)
        .flat()
        .filter((checkbox) => checkbox !== undefined);
    const uncheckAll = () => {
        headings.forEach((heading) => {
            if (hCheckboxes) {
                hCheckboxes.forEach((checkbox) => {
                    const checkboxElement = document.getElementById(checkbox.id);
                    checkboxElement.checked = false;
                    Checkboxes.set([]);
                });
                hCheckboxes = hCheckboxes; // trigger reactivity
            }
        });
    };
    onMount(() => {
        // set the first heading as active on initial load
        if (!$currentHeading || !headings.find((h) => h.slug === $currentHeading)) {
            currentHeading.set(headings[0]?.slug);
        }
        currentHeading.subscribe((slug) => {
            //check if any heading has the same slug as the currentHeading
            const heading = headings.find((h) => h.slug === slug);
            activeHeading = heading?.slug || $currentHeading;
            // wait 1 second for sidebar selection animation to finish

            const active = document.querySelector(".toc nav-item.active");
            if (active) {
                active.scrollIntoView({ block: "nearest" });
            }
        });
        if (hCheckboxes) {
            Checkboxes.subscribe((checks) => {
                hCheckboxes.forEach((hCheckbox) => {
                    hCheckbox.checked = checks.find((check) => check.id === hCheckbox.id)?.checked || false;
                });
                hCheckboxes = hCheckboxes; // trigger reactivity
            });
        }
    });
</script>

<div class="nav flex-column sticky-top-under align-items-end pt-1">
    <div class="d-none d-xl-block w-100">
        {@render rightSidebarTop()}
        {#if showToc}
            <strong class="h5 my-2 text-body">On this page</strong>
        {/if}
        <nav id="TableOfContents" class="d-flex flex-column">
            {#if showToc}
                <ul class="mb-0 mt-1">
                    {#each headings as heading (heading)}
                        <li
                            class={"nav-item ms-1 heading-padding-" + (heading.depth - minHeadingDepth)}
                            class:active={heading.slug === activeHeading}
                            class:collapse={heading.hidden && !$showHidden}
                        >
                            <a
                                class="nav-link py-1 ps-3 text-body-secondary small d-inline-flex align-items-center"
                                href={"#" + heading.slug}
                            >
                                {#if heading.fa_icon}
                                    <i class={heading.fa_icon + " fa-fw me-2"} aria-hidden="true"></i>
                                {/if}
                                {@html heading.text}
                                {#if hCheckboxes.find((hc) => hc?.id.startsWith("checkbox-" + heading.slug))}
                                    <span class="ms-2">
                                        <ProgressIndicator
                                            progress={(hCheckboxes.filter(
                                                (check) =>
                                                    check?.id.startsWith("checkbox-" + heading.slug) && check?.checked,
                                            ).length /
                                                hCheckboxes.filter((check) =>
                                                    check?.id.startsWith("checkbox-" + heading.slug),
                                                ).length) *
                                                100}
                                            size={25}
                                            strokeWidth={4}
                                            isCurrent={true}
                                            confetti={true}
                                        />
                                    </span>
                                {/if}
                            </a>
                        </li>
                    {/each}
                </ul>

                {#if hCheckboxes.length > 0}
                    <button class="btn btn-sm btn-outline-secondary mt-2" on:click={uncheckAll} on:keydown={uncheckAll}>
                        Uncheck all
                    </button>
                {/if}

                <div class="d-flex flex-column flex-xl-row justify-content-between mt-xl-2">
                    <a
                        href="#/"
                        class="back-to-top text-body-secondary text-small mb-2"
                        on:click={() => window.scrollTo(0, 0)}
                    >
                        <i class="fa-solid fa-arrow-up-to-line" aria-hidden="true"></i> Back to top
                    </a>
                    {@render rightSidebarLinkBar()}
                </div>
            {/if}
            {@render defaultContent()}
        </nav>
    </div>
</div>

<style lang="scss">
    nav > ul {
        font-size: 0.875rem;
        list-style: none;
        padding-left: 0;
    }
    .nav ul {
        overflow-y: auto;
        max-height: calc(100vh - 8rem);
    }

    li {
        border-inline-start: 2pt solid var(--bs-border-color);
        transition:
            background-color 0.3s ease-out,
            border-left 0.3s ease-out;
        scroll-margin-top: 6rem;
        scroll-margin-bottom: 6rem;
        &:hover {
            background-color: rgba(var(--bs-success), 0.85);
        }
    }

    li.active {
        border-left: 2pt solid var(--bs-success);
        background-color: rgba(var(--bs-success), 0.85);
        a {
            color: var(--bs-green-800) !important;
        }
    }

    :global([data-bs-theme="dark"]) {
        li {
            &:hover {
                & a {
                    color: var(--bs-gray-300) !important;
                }
            }
        }
        li.active {
            & a {
                color: var(--bs-gray-300) !important;
            }
        }
    }
    @mixin heading-padding($number) {
        @if $number >= 0 and $number <= 6 {
            padding-left: #{$number * 1rem};
        } @else {
            @warn "Invalid input. Please provide a number between 1 and 6.";
        }
    }

    @for $i from 0 through 6 {
        .heading-padding-#{$i} {
            @include heading-padding($i);
        }
    }
</style>
