<script lang="ts">
    import docsearch from "@docsearch/js";

    import "@docsearch/css";
    import { onMount } from "svelte";

    interface Props {
        id?: string;
    }

    let { id = "" }: Props = $props();
    onMount(() => {
        docsearch({
            container: "#" + id,
            appId: "01BY3A8NRJ",
            indexName: "nf-co",
            apiKey: "c726615ab69f88b4e26bc891c97c1808",
            placeholder: "Search",
        });
    });
</script>

<div {id}></div>

<style lang="scss">
    :global(:root) {
        --docsearch-text-color: var(--bs-body-color);
        --docsearch-searchbox-background: var(--bs-body-bg);
        --docsearch-searchbox-shadow: none;
        --docsearch-searchbox-focus-background: var(--bs-body-bg);
        --docsearch-primary-color: #{var(--bs-success)};
        --docsearch-container-background: #{rgba(var(--bs-black), 0.5)};
        --docsearch-modal-background: #{var(--bs-body-bg)};
        --docsearch-modal-shadow: "none";
        --docsearch-footer-shadow: none;
        --docsearch-muted-color: #{var(--bs-secondary)};
        --docsearch-secondary-text-color: #{var(--bs-tertiary-emphasis)};
        --docsearch-icon-color: #{var(--bs-tertiary-emphasis)};
        --docsearch-highlight-color: #{var(--bs-tertiary-emphasis)};
    }
    :global([data-bs-theme="dark"]:root) {
        --docsearch-hit-color: #{var(--bs-body-color)};
        --docsearch-icon-color: #{var(--bs-body-color)};
        --docsearch-logo-color: #{var(--bs-body-color)};
        --docsearch-hit-background: #{var(--bs-body-bg)};
        --docsearch-hit-shadow: "none";
        --docsearch-footer-background: #{var(--bs-tertiary-bg)};
        --docsearch-muted-color: #{var(--bs-tertiary-emphasis)};
        --docsearch-subtle-color: #{var(--bs-tertiary-emphasis)};
    }
    :global(.DocSearch-Form) {
        border: 1pt solid var(--bs-border-color);
    }
    :global(button.DocSearch-Button) {
        border-radius: var(--bs-border-radius);
        margin-left: 0;
    }
    :global(button.DocSearch-Button .DocSearch-Button-Keys) {
        background: none;
        border-radius: var(--bs-border-radius-sm);

        //   padding: 0.2rem;
        margin-top: 0.15rem;

        min-width: 3rem;
        justify-content: space-around;
        align-items: center;

        display: none;
        /*
        @include media-breakpoint-up(xxl) {
            display: flex;
        }
        */
    }
    :global(.DocSearch-Button-Key, .DocSearch-Commands-Key) {
        align-items: center;
        background: var(--bs-body-bg);
        box-shadow: 0 0 0 1px var(--bs-secondary);
        margin-right: 0.2rem;
        padding: 0;
    }
    :global(.DocSearch .DocSearch-Hits) {
        &:last-of-type {
            margin-bottom: 0;
        }
        &::before {
            height: 0;
        }
    }
    :global(.DocSearch-Hit) {
        & :global(a) {
            box-shadow: none;
        }
        & :global(mark) {
            padding: 0;
        }
    }
    :global(.DocSearch-Hit-icon:has(svg[height="20"])) {
        line-height: 20px;
    }

    :global(.DocSearch-Hit-source) {
        font-size: 1rem;
    }
    :global([data-bs-theme="dark"] .DocSearch-Hit[aria-selected="true"] a) {
        background-color: rgba(var(--bs-success-rgb), 0.75);
    }
    :global([data-bs-theme="dark"]) {
        :global(.cls-1),
        :global(.cls-2) {
            fill: #{var(--bs-body-color)};
        }
    }
</style>
