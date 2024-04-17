<script lang="ts">
    import { currentHeading, showHidden, checkboxes } from '@components/store';
    import { onMount } from 'svelte';
    import ProgressIndicator from '@components/sidebar/ProgressIndicator.svelte';

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

            const active = document.querySelector('.toc nav-item.active');
            if (active) {
                active.scrollIntoView({ block: 'nearest' });
            }
        });
    });

    checkboxes.subscribe((checkboxes) => {
        checkboxes.forEach((checkbox) => {
            headings.forEach((heading) => {
                if (heading.checkboxes) {
                    heading.checkboxes.forEach((hCheckbox) => {
                        if (hCheckbox.id === checkbox.id) {
                            hCheckbox.checked = checkbox.checked;
                        }
                    });
                }
            });
        });
    });

    const uncheckAll = () => {
        headings.forEach((heading) => {
            if (heading.checkboxes) {
                heading.checkboxes.forEach((checkbox) => {
                    const checkboxElement = document.getElementById(checkbox.id);
                    checkboxElement.checked = false;
                    checkboxes.set([]);
                });
            }
        });
    };
</script>

<div class="nav flex-column sticky-top-under align-items-end pt-1">
    <div class="d-none d-md-block w-100">
        <slot name="right-sidebar-top" />
        {#if showToc}
            <strong class="h6 my-2 text-body">On this page</strong>
        {/if}
        <!-- <hr class="my-1" /> -->
        <nav id="TableOfContents" class="d-none d-md-flex flex-column">
            {#if showToc}
                <ul class="mb-0 mt-1">
                    {#each headings as heading (heading)}
                        <li
                            class={'nav-item heading-padding-' + (heading.depth - minHeadingDepth)}
                            class:active={heading.slug === activeHeading}
                            class:collapse={heading.hidden && !$showHidden}
                        >
                            <a class="nav-link py-1 ps-3 text-body-secondary small" href={'#' + heading.slug}>
                                {#if heading.fa_icon}
                                    <i class={heading.fa_icon + ' fa-fw '} aria-hidden="true" />
                                {/if}
                                {@html heading.text}
                                {#if heading.checkboxes}
                                    <ProgressIndicator
                                        progress={(heading.checkboxes.filter((check) => check.checked).length /
                                            heading.checkboxes.length) *
                                            100}
                                    />
                                {/if}
                            </a>
                        </li>
                    {/each}
                </ul>
                <div class="">
                    <a
                        href="#/"
                        class="back-to-top text-body-secondary text-small float-end mb-2"
                        on:click={() => window.scrollTo(0, 0)}
                    >
                        <i class="fa-solid fa-arrow-up-to-line" aria-hidden="true" /> Back to top
                    </a>
                </div>
            {/if}
            {#if headings.some((h) => h.checkboxes)}
                <button class="btn btn-sm btn-outline-secondary" on:click={uncheckAll} on:keydown={uncheckAll}>
                    Uncheck all
                </button>
            {/if}
            <slot />
        </nav>
    </div>
</div>

<style lang="scss">
    @import 'src/styles/_variables.scss';

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
        border-inline-start: 2pt solid $border-color;
        transition:
            background-color 0.3s ease-out,
            border-left 0.3s ease-out;
        scroll-margin-top: 6rem;
        scroll-margin-bottom: 6rem;
        &:hover {
            background-color: transparentize($success, 0.85);
        }
    }

    li.active {
        a {
            color: $green-800 !important;
        }
        border-left: 2pt solid $success;
        background-color: transparentize($success, 0.75);
    }

    :global([data-bs-theme='dark']) {
        li {
            border-inline-start: 2pt solid $border-color-dark;
            &:hover {
                background-color: transparentize($success-dark, 0.6);
                & a {
                    color: $gray-300 !important;
                }
            }
        }
        li.active {
            border-left: 2pt solid $success-dark;
            background-color: transparentize($success-dark, 0.35);
            & a {
                color: $gray-300 !important;
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
