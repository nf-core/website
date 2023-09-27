<script lang="ts">
    import docsearch from '@docsearch/js';

    import '@docsearch/css';
    import { onMount } from 'svelte';

    export let id: string = '';
    onMount(() => {
        docsearch({
            container: '#' + id,
            appId: '01BY3A8NRJ',
            indexName: 'nf-co',
            apiKey: 'c726615ab69f88b4e26bc891c97c1808',
            placeholder: 'Search',
        });
    });
</script>

<div {id} />

<style lang="scss">
    @import '@styles/_variables.scss';

    :global(:root) {
        --docsearch-text-color: #{$input-color};
        --docsearch-searchbox-background: #{$input-bg};
        --docsearch-searchbox-shadow: none;
        --docsearch-searchbox-focus-background: #{$input-bg};
        --docsearch-primary-color: #{$nf-core-color};
        --docsearch-container-background: #{rgba($black, 0.5)};
        --docsearch-modal-shadow: 'none';
        --docsearch-footer-shadow: none;
    }
    :global(.DocSearch-Modal) {
        background: $body-secondary-bg;
    }
    :global(button.DocSearch-Button) {
        border-radius: $border-radius;
        margin-left: 0;
    }
    :global(button.DocSearch-Button .DocSearch-Button-Keys) {
        background: none;
        border-radius: $border-radius-sm;

        //   padding: 0.2rem;
        margin-top: 0.15rem;

        min-width: 3rem;
        justify-content: space-around;
        align-items: center;
    }
    :global(.DocSearch-Button-Key, .DocSearch-Commands-Key) {
        align-items: center;
        background: $white;
        box-shadow: 0 0 0 1px $secondary;
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
    :global(.DocSearch-Hit-icon:has(svg[height='20'])) {
        line-height: 20px;
    }

    :global([data-bs-theme='dark']:root) {
        --docsearch-text-color: #{$input-color};
        --docsearch-searchbox-background: #{$input-bg};
        --docsearch-searchbox-shadow: none;
        --docsearch-searchbox-focus-background: #{$input-bg};
        --docsearch-primary-color: #{$nf-core-color-dark};
        --docsearch-container-background: #{rgba($black, 0.5)};
        --docsearch-modal-shadow: 'none';
        --docsearch-hit-color: #{$body-color-dark};
        --docsearch-hit-background: #{$body-bg-dark};
        --docsearch-hit-shadow: 'none';
        --docsearch-footer-background: #{$body-secondary-bg-dark};
        --docsearch-footer-shadow: none;
    }
    :global(.DocSearch-Hit-source) {
        font-size: 1rem;
    }
    :global([data-bs-theme='dark'] .DocSearch-Modal, [data-bs-theme='dark'] .DocSearch-Hit-source) {
        background: $body-secondary-bg-dark;
    }
    :global(.DocSearch-Container) {
        @include media-breakpoint-down(lg) {
            padding-top: 3rem;
        }
    }
    :global(.double-navbar ~ .DocSearch-Container) {
        @include media-breakpoint-down(lg) {
            padding-top: 5rem;
        }
    }
    :global([data-bs-theme='dark'] .DocSearch-Button-Key, [data-bs-theme='dark'] .DocSearch-Commands-Key) {
        background: $black;
        box-shadow: 0 0 0 1px $secondary-border-subtle-dark;
    }
</style>
