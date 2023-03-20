<script lang="ts">
    const icon_mdi_dots_vertical = `<svg viewBox="0 0 24 24" class="d-inline-block" astro-icon="mdi:dots-vertical"><path fill="currentColor" d="M12 16a2 2 0 0 1 2 2 2 2 0 0 1-2 2 2 2 0 0 1-2-2 2 2 0 0 1 2-2m0-6a2 2 0 0 1 2 2 2 2 0 0 1-2 2 2 2 0 0 1-2-2 2 2 0 0 1 2-2m0-6a2 2 0 0 1 2 2 2 2 0 0 1-2 2 2 2 0 0 1-2-2 2 2 0 0 1 2-2z"></path></svg>`;
    const icon_mdi_close = `<svg viewBox="0 0 24 24" class="d-inline-block" astro-icon="mdi:close"><path fill="currentColor" d="M19 6.41 17.59 5 12 10.59 6.41 5 5 6.41 10.59 12 5 17.59 6.41 19 12 13.41 17.59 19 19 17.59 13.41 12 19 6.41z"></path></svg>`;

    let visible = false;
    function toggleVisible() {
        visible = !visible;
    }
</script>

<header class="docs-nav small border-bottom d-md-none">
    <div class="me-2 w-100 overflow-x-hidden text-truncate text-nowrap">
        <button class="btn" on:click={toggleVisible}>
            {@html icon_mdi_dots_vertical}
        </button>
        <slot name="title" />
    </div>
</header>

{#if visible}
    <div class="fixed inset-0 z-50">
        <div class="fixed inset-0 backdrop-blur-sm" />
        <nav class="side-nav relative h-screen w-80 max-w-full overflow-y-auto bg-gray-800 p-5 text-gray-400">
            <button class="float-right mt-1" on:click={toggleVisible}>
                {@html icon_mdi_close}
            </button>
            <h4 class="mb-4 text-base font-semibold">MultiQC Documentation</h4>
            <slot name="menu" />
        </nav>
    </div>
{/if}

<style lang="scss">
    header.small {
        margin-top: 3rem;
    }
    :global(body:has(.docs-nav) .content) {
        padding-top: 0;
    }
    button svg {
        width: 1.5rem;
        height: 1.5rem;
    }
</style>
