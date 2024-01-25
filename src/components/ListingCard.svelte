<script lang="ts">
    export let recentRelease: boolean = false;
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
</div>

<style lang="scss">
    .card {
        width: 100%;
        height: 100%;
    }
</style>
