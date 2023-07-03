<script lang="ts">
    import { currentPage } from '@components/store';
    import { onMount } from 'svelte';

    export let lastPage: number = 1;

    let pages = [];
    let truncatedPages = [];

    const maxPages = 7;
    let truncated = false;

    onMount(() => {
        generatePages();
    });

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
            <li
                class="page-item"
                on:click={() => handlePageChange($currentPage - 1)}
                on:keydown={() => handlePageChange($currentPage - 1)}
                class:disabled={$currentPage === 1}
                aria-disabled={$currentPage === 1}
            >
                <a class="page-link" href="#" tabindex="-1">Previous</a>
            </li>

            {#if truncated}
                {#if truncatedPages[0] > 1}
                    <li class="page-item">
                        <a class="page-link" href="#">1</a>
                    </li>
                    <li class="page-item disabled">
                        <a class="page-link" href="#">...</a>
                    </li>
                {/if}
                {#each truncatedPages as page}
                    <li
                        class="page-item"
                        on:click={() => handlePageChange(page)}
                        on:keydown={() => handlePageChange(page)}
                        class:active={$currentPage === page}
                    >
                        <a class="page-link" href="#">{page}</a>
                    </li>
                {/each}
                {#if truncatedPages[truncatedPages.length - 1] < lastPage - 1}
                    <li class="page-item disabled">
                        <a class="page-link" href="#">...</a>
                    </li>
                    <li
                        class="page-item"
                        on:click={() => handlePageChange(lastPage)}
                        on:keydown={() => handlePageChange(lastPage)}
                        class:active={$currentPage === lastPage}
                    >
                        <a class="page-link" href="#">{lastPage}</a>
                    </li>
                {/if}
            {:else}
                {#each pages as page}
                    <li
                        class="page-item"
                        on:click={() => handlePageChange(page)}
                        on:keydown={() => handlePageChange(page)}
                        class:active={$currentPage === page}
                    >
                        <a class="page-link" href="#">{page}</a>
                    </li>
                {/each}
            {/if}

            <li
                class="page-item"
                on:click={() => handlePageChange($currentPage + 1)}
                on:keydown={() => handlePageChange($currentPage + 1)}
                class:disabled={$currentPage === lastPage}
                aria-disabled={$currentPage === lastPage}
            >
                <a class="page-link" href="#">Next</a>
            </li>
        </ul>
    </nav>
</div>
