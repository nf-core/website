<script lang="ts">
    let visible = true;
    function toggleVisible() {
        visible = !visible;
    }
</script>

<div class="docs-nav fixed-top bg-body small border-bottom d-md-none">
    <div class="w-100 text-nowrap">
        <button class="btn text-body d-flex align-items-center ps-2" on:click={toggleVisible}>
            <i class="fa-regular fa-ellipsis-vertical me-2" />

            <slot name="title" />
        </button>
    </div>
</div>

{#if visible}
    <div class="d-md-none position-fixed bg-body min-vh-100 z-3">
        <span
            class="position-fixed bg-dark bg-opacity-50 w-100 min-vh-100"
            on:click={toggleVisible}
            on:keypress={toggleVisible}
        />
        <nav class="side-nav position-relative bg-body py-3 px-0 text-gray-400 overflow-y-auto">
            <button type="button" class="btn-close float-end me-2" on:click={toggleVisible} aria-label="Close" />
            <slot name="menu" />
        </nav>
    </div>
{/if}

<style lang="scss">
    @import '@styles/_variables.scss';
    .docs-nav {
        margin-top: 3rem;
        z-index: 1029; // reduce z-index by one to have navbarToc on top
    }
    .side-nav {
        margin-top: 5rem;
        max-height: calc(100vh - 5rem);
        width: 75vw;
    }
</style>
