<script lang="ts">
    export let recentRelease: boolean;
</script>

<div class="card m-2">
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
    @import '@styles/_variables.scss';
    .card {
        width: 32rem;
    }
    @include media-breakpoint-down(md) {
        .card {
            width: 100%;
        }
    }
</style>
