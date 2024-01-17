<script lang="ts">
    import { currentPage } from '@components/store';
    import { onMount } from 'svelte';

    export let lastPage: number = 1;

    let pages: number[] = [];
    let truncatedPages: number[] = [];

    const maxPages = 7;
    let truncated = false;

    $: {
        if (lastPage > 0) {
            generatePages();
        }
    }

    function generatePages() {
        pages = [];
        truncatedPages = [];
        truncated = lastPage > maxPages;

        if (truncated) {
            const startIndex = Math.max($currentPage - Math.floor(maxPages / 2), 1);
            const endIndex = Math.min(startIndex + maxPages - 1, lastPage);

            for (let i = startIndex; i <= endIndex; i++) {
                truncatedPages.push(i);
            }
        } else {
            for (let i = 1; i <= lastPage; i++) {
                pages.push(i);
            }
        }
    }

    function handlePageChange(page) {
        $currentPage = page;
        generatePages();
    }
</script>

<div class="d-flex justify-content-center mt-2">
    <nav aria-label="Page navigation">
        <ul class="pagination">
            <li class="page-item" class:disabled={$currentPage === 1}>
                <span
                    on:click={() => handlePageChange($currentPage - 1)}
                    on:keydown={() => handlePageChange($currentPage - 1)}
                    role="button"
                    class="page-link"
                    tabindex="0">Previous</span
                >
            </li>

            {#if truncated}
                {#if truncatedPages[0] > 1}
                    <li class="page-item">
                        <span class="page-link">1</span>
                    </li>
                    <li class="page-item disabled">
                        <span class="page-link">...</span>
                    </li>
                {/if}
                {#each truncatedPages as page}
                    <li class="page-item" class:active={$currentPage === page}>
                        <span
                            on:click={() => handlePageChange(page)}
                            on:keydown={() => handlePageChange(page)}
                            role="button"
                            class="page-link"
                            tabindex="0">{page}</span
                        >
                    </li>
                {/each}
                {#if truncatedPages[truncatedPages.length - 1] < lastPage - 1}
                    <li class="page-item disabled">
                        <span class="page-link">...</span>
                    </li>
                    <li class="page-item" class:active={$currentPage === lastPage}>
                        <span
                            on:click={() => handlePageChange(lastPage)}
                            on:keydown={() => handlePageChange(lastPage)}
                            role="button"
                            tabindex="0"
                            class="page-link">{lastPage}</span
                        >
                    </li>
                {/if}
            {:else}
                {#each pages as page}
                    <li class="page-item" class:active={$currentPage === page}>
                        <span
                            on:click={() => handlePageChange(page)}
                            on:keydown={() => handlePageChange(page)}
                            role="button"
                            tabindex="0"
                            class="page-link">{page}</span
                        >
                    </li>
                {/each}
            {/if}

            <li class="page-item" class:disabled={$currentPage === lastPage}>
                <span
                    on:click={() => handlePageChange($currentPage + 1)}
                    on:keydown={() => handlePageChange($currentPage + 1)}
                    class="page-link"
                    role="button"
                    tabindex="0">Next</span
                >
            </li>
        </ul>
    </nav>
</div>
