<script>
    import { EventIsOngoing } from '@components/store';
</script>

<div class="event-indicator align-center d-inline-block" class:d-none={!$EventIsOngoing}>
    <slot />
</div>

<style lang="scss">
    .event-indicator :global(svg[data-icon]) {
        color: var(--bs-danger);
        width: 10pt;
        margin-bottom: 0.1rem;
    }
    .event-indicator :global(svg[data-icon] .fa-secondary) {
        opacity: 0.4;
    }
</style>
