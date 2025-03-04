<script lang="ts">
    interface Props {
        title?: import("svelte").Snippet;
        menu?: import("svelte").Snippet;
    }

    let { title, menu }: Props = $props();
    let visible = $state(false);
    function toggleVisible() {
        visible = !visible;
    }
</script>

<div class="docs-nav fixed-top bg-body small border-bottom d-md-none">
    <div class="w-100 text-nowrap">
        <button class="btn text-body d-flex align-items-center ps-2" onclick={toggleVisible}>
            <i class="fa-regular fa-ellipsis-vertical me-2"></i>

            {@render title?.()}
        </button>
    </div>
</div>

{#if visible}
    <div class="d-md-none bg-body z-3">
        <span
            class="position-fixed bg-dark bg-opacity-75 w-100 min-vh-100"
            onclick={toggleVisible}
            onkeypress={toggleVisible}
            role="button"
            tabindex="0"
        ></span>
        <nav class="side-nav bg-body pb-2 px-0 text-gray-400 overflow-y-auto">
            {@render menu?.()}
        </nav>
    </div>
{/if}

<style lang="scss">
    .docs-nav {
        margin-top: 3rem;
    }
    .side-nav {
        padding-top: 3rem; //offset for navbar + sectionNavbar from sidebarNav
        max-height: calc(100dvh - 3rem);
        width: 100vw;
    }
    :global(.navbar.fixed-top:has(.collapse.show)) {
        & .docs-nav,
        & .side-nav {
            display: none;
        }
    }
</style>
