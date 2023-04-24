<script lang="ts">
    const icon_mdi_dots_vertical = `<svg viewBox="0 0 24 24" class="d-inline-block" astro-icon="mdi:dots-vertical"><path fill="currentColor" d="M12 16a2 2 0 0 1 2 2 2 2 0 0 1-2 2 2 2 0 0 1-2-2 2 2 0 0 1 2-2m0-6a2 2 0 0 1 2 2 2 2 0 0 1-2 2 2 2 0 0 1-2-2 2 2 0 0 1 2-2m0-6a2 2 0 0 1 2 2 2 2 0 0 1-2 2 2 2 0 0 1-2-2 2 2 0 0 1 2-2z"></path></svg>`;

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
    <div class="d-md-none position-fixed bg-body min-vh-100 z-3">
        <span
            class="position-fixed bg-dark bg-opacity-50 w-100 min-vh-100"
            on:click={toggleVisible}
            on:keypress={toggleVisible}
        />
        <nav class="side-nav position-relative w-100 bg-body p-3 pe-0 text-gray-400 overflow-y-auto">
            <button type="button" class="btn-close float-end me-2" on:click={toggleVisible} aria-label="Close" />
            <h4 class="mb-2 fw-semibold">nf-core Documentation</h4>
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
    @include media-breakpoint-down(md) {
        :global(body:has(.docs-nav) .content) {
            padding-top: 5.5rem;
        }
    }
    .side-nav {
        margin-top: 5rem;
        max-height: calc(100vh - 5rem);
    }
</style>
