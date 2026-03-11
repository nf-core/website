<script lang="ts">
    interface Props {
        type: "event" | "blog";
        timeSpan?: [number, number];
        children?: import("svelte").Snippet;
    }

    let { type, timeSpan = [0, 0], children }: Props = $props();

    const now = new Date().getTime();
    const show = $derived(timeSpan[0] < now && timeSpan[1] > now);
</script>

<div class={`event-indicator align-center d-inline-block ${type}`} class:d-none={!show}>
    {@render children?.()}
</div>

<style lang="scss">
    .event-indicator :global(svg[data-icon]) {
        width: 10pt;
        margin-bottom: 0.1rem;
    }
    .event-indicator :global(svg[data-icon] .fa-secondary) {
        opacity: 0.4;
    }
    .event-indicator.event {
        color: var(--bs-danger);
    }
    .event-indicator.blog {
        color: var(--bs-primary);
    }
</style>
