<script lang="ts">
    export let recentRelease: boolean = false;
    export let footer: boolean = false;
    export let cardHeader: () => any = () => null;
    export let cardBody: () => any = () => null;
    export let cardFooter: () => any = () => null;
</script>

<div class="card">
    <div class="card-header" class:border-success={recentRelease}>
        <h2 class="mb-0">
            {@render cardHeader()}
        </h2>
    </div>
    <div class="card-body" class:pt-0={recentRelease}>
        {@render cardBody()}
    </div>
    {#if footer}
        <div class="card-footer h-100">
            {@render cardFooter()}
        </div>
    {/if}
</div>

<style lang="scss">
    .card {
        width: 100%;
        height: 100%;
    }
</style>
