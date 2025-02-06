<script lang="ts">
    import { run } from 'svelte/legacy';

    import { currentPage } from '@components/store';
    import { onMount } from 'svelte';

    interface Props {
        lastPage?: number;
    }

    let { lastPage = 1 }: Props = $props();

    let pages: number[] = $state([]);
    let truncatedPages: number[] = $state([]);

    const maxPages = 7;
    let truncated = $state(false);

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
    run(() => {
        if (lastPage > 0) {
            generatePages();
        }
    });
</script>

<div class="d-flex justify-content-center mt-2">
    <nav aria-label="Page navigation">
        <ul class="pagination">
            <li class="page-item" class:disabled={$currentPage === 1}>
                <span
                    onclick={() => handlePageChange($currentPage - 1)}
                    onkeydown={() => handlePageChange($currentPage - 1)}
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
                            onclick={() => handlePageChange(page)}
                            onkeydown={() => handlePageChange(page)}
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
                            onclick={() => handlePageChange(lastPage)}
                            onkeydown={() => handlePageChange(lastPage)}
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
                            onclick={() => handlePageChange(page)}
                            onkeydown={() => handlePageChange(page)}
                            role="button"
                            tabindex="0"
                            class="page-link">{page}</span
                        >
                    </li>
                {/each}
            {/if}

            <li class="page-item" class:disabled={$currentPage === lastPage}>
                <span
                    onclick={() => handlePageChange($currentPage + 1)}
                    onkeydown={() => handlePageChange($currentPage + 1)}
                    class="page-link"
                    role="button"
                    tabindex="0">Next</span
                >
            </li>
        </ul>
    </nav>
</div>
