<script lang="ts">
    export let type: 'event' | 'blog';
    export let timeSpan: [number, number] = [0, 0];

    const now = new Date().getTime();
    const show = timeSpan[0] < now && timeSpan[1] > now;
</script>

<div class={`event-indicator align-center d-inline-block ${type}`} class:d-none={!show}>
    <slot />
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
