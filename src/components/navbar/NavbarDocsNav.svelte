<script lang="ts">
    let visible = false;
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
    <div class="d-md-none bg-body z-3">
        <span
            class="position-fixed bg-dark bg-opacity-75 w-100 min-vh-100"
            on:click={toggleVisible}
            on:keypress={toggleVisible}
            role="button"
            tabindex="0"
        />
        <nav class="side-nav bg-body pb-2 px-0 text-gray-400 overflow-y-auto">
            <button type="button" class="btn-close float-end me-2" on:click={toggleVisible} aria-label="Close" />
            <slot name="menu" />
        </nav>
    </div>
{/if}

<style lang="scss">
    @import '../../styles/_variables.scss';
    .docs-nav {
        margin-top: 3rem;
    }
    .side-nav {
        padding-top: 3rem; //offset for navbar + sectionNavbar from sidebarNav
        max-height: calc(100dvh - 3rem);
        width: 75vw;
    }
    :global(.navbar.fixed-top:has(.collapse.show)) .docs-nav {
        display: none;
    }
</style>
