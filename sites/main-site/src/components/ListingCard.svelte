<script lang="ts">
    import type { Snippet } from "svelte";

    interface Props {
        recentRelease?: boolean;
        footer?: boolean;
        cardClass?: string;
        cardStyle?: string;
        cardHeader?: Snippet;
        cardBody?: Snippet;
        cardFooter?: Snippet;
    }

    let {
        recentRelease = false,
        footer = false,
        cardClass = "",
        cardStyle = "",
        cardHeader,
        cardBody,
        cardFooter,
    }: Props = $props();
</script>

<div class={["card", cardClass]} style={cardStyle}>
    <div class="card-header" class:border-success={recentRelease}>
        <h2 class="mb-0">
            {@render cardHeader?.()}
        </h2>
    </div>
    <div class="card-body" class:pt-0={recentRelease}>
        {@render cardBody?.()}
    </div>
    {#if footer}
        <div class="card-footer h-100">
            {@render cardFooter?.()}
        </div>
    {/if}
</div>

<style lang="scss">
    .card {
        width: 100%;
        &:global(.rounded-start-0) {
            border-left: 5px solid;
        }
    }
</style>
