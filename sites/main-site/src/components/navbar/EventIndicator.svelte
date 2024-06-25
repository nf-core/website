<script lang="ts">
    import { EventIsOngoing, newBlogPost } from '@components/store';

    export let type: 'event' | 'blog';

    let show = false;

    if (type === 'event' && EventIsOngoing) {
        show = true;
    } else if (type === 'blog' && newBlogPost) {
        show = true;
    }
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
