<script>
    import { EventIsOngoing } from '@components/store';
</script>

{#if $EventIsOngoing}
    <div class="event-indicator align-center">
        <slot />
    </div>
{/if}

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
