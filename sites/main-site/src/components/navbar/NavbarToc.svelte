<script lang="ts">
    import { currentHeading } from "@components/store";
    import SchemaListingTocButtons from "@components/schema/SchemaListingTocButtons.svelte";
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
</script>

<div class="d-md-none toc-md">
    <div class="dropdown">
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
                    <a class="dropdown-item" href={"#" + heading.slug}>
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
            overflow-y: auto;
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
