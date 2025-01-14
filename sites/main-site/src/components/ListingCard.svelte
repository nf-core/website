<!-- @migration-task Error while migrating Svelte code: This migration would change the name of a slot making the component unusable -->
<script lang="ts">
    export let recentRelease: boolean = false;
    export let footer: boolean = false;
</script>

<div class="card">
    <div class="card-header" class:border-success={recentRelease}>
        <h2 class="mb-0">
            <slot name="card-header" />
        </h2>
    </div>
    <div class="card-body" class:pt-0={recentRelease}>
        <slot name="card-body" />
    </div>
    {#if footer}
        <div class="card-footer h-100">
            <slot name="card-footer" />
        </div>
    {/if}
</div>

<style lang="scss">
    .card {
        width: 100%;
        height: 100%;
    }
</style>
