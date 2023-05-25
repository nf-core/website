<script>
    import { EventIsOngoing } from '@components/store';
</script>

{#if $EventIsOngoing}
    <span class="fa-stack small text-danger">
        <i class="fa-duotone fa-circle fa-stack-1x" />
        <i class="fas fa-circle-small fa-stack-1x" />
    </span>
{/if}
