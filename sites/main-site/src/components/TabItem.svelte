<script lang="ts">
    import { CurrentTab } from "@components/store";
    interface Props {
        name: string;
        children?: import("svelte").Snippet;
    }

    let { name, children }: Props = $props();
</script>

<div
    class={`tab-pane fade border border-top-0 p-2 ${name}-tab-pane`}
    id={`${name}-tab-pane`}
    class:active={name === $CurrentTab}
    class:show={name === $CurrentTab}
    class:collapse={name !== $CurrentTab}
    role="tabpanel"
    aria-labelledby={`${name}-tab`}
    tabindex="0"
>
    {@render children?.()}
</div>

<style lang="scss">
    .tab-pane {
        margin-bottom: 1rem;
        :global(.admonition),
        :global(figure),
        :global(pre) {
            &:last-child {
                margin-bottom: 0.25rem !important;
            }
        }
    }
</style>
