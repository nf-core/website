<script>
    import { EventIsOngoing } from '@components/store';
</script>

{#if $EventIsOngoing}
    <div class="event-indicator align-center">
        <slot />
    </div>
{/if}

<style lang="scss">
    @import '@styles/_variables.scss';
    .event-indicator :global(svg[astro-icon]) {
        fill: $danger;
        width: 10pt;
        margin-bottom: 0.1rem;
    }
    .event-indicator :global(svg[astro-icon] .fa-secondary) {
        opacity: 0.4;
    }
</style>
