<script lang="ts">
    import PipelineCard from '@components/pipeline/PipelineCard.svelte';
    import { CurrentFilter, Filters, SortBy, DisplayStyle, SearchQuery } from '@components/store';

    export let pipelines: {
        name: string;
        description: string;
        stargazers_count: number;
        topics: string[];
        releases: {
            published_at: string;
            tag_name: string;
        }[];
        archived: boolean;
    }[] = [];

    let sortInverse = false;
    function handleSort(sor) {
        if (sor === $SortBy) {
            sortInverse = !sortInverse;
        } else {
            sortInverse = false;
        }
        SortBy.set(sortInverse ? sor + ';inverse' : sor);
    }
    const searchPipelines = (pipeline) => {
        if ($SearchQuery === '') {
            return true;
        }
        if (pipeline.name.toLowerCase().includes($SearchQuery.toLowerCase())) {
            return true;
        }
        if (pipeline.description && pipeline.description.toLowerCase().includes($SearchQuery.toLowerCase())) {
            return true;
        }
        if (pipeline.topics.some((topic) => topic.toLowerCase().includes($SearchQuery.toLowerCase()))) {
            return true;
        }
        return false;
    };

    const filterPipelines = (pipeline) => {
        if ($CurrentFilter.includes('Released') && pipeline.releases.length > 1 && !pipeline.archived) {
            return true;
        }
        if ($CurrentFilter.includes('Under development') && pipeline.releases.length === 1 && !pipeline.archived) {
            return true;
        }
        if ($CurrentFilter.includes('Archived') && pipeline.archived === true) {
            return true;
        }
        return false;
    };

    const sortPipelines = (a, b) => {
        sortInverse = $SortBy.endsWith(';inverse');
        if ($SortBy.startsWith('Alphabetical')) {
            if (sortInverse) {
                return b.name.localeCompare(a.name);
            } else {
                return a.name.localeCompare(b.name);
            }
        } else if ($SortBy === 'Stars') {
            if (sortInverse) {
                return a.stargazers_count - b.stargazers_count;
            } else {
                return b.stargazers_count - a.stargazers_count;
            }
        } else if ($SortBy === 'Last release') {
            // handle case where a pipeline has no releases
            if (a.releases.length === 1) {
                return 1 * (sortInverse ? -1 : 1);
            }
            if (b.releases.length === 1) {
                return -1 * (sortInverse ? -1 : 1);
            }
            if (sortInverse) {
                return new Date(a.releases[0].published_at) - new Date(b.releases[0].published_at);
            } else {
                return new Date(b.releases[0].published_at) - new Date(a.releases[0].published_at);
            }
        }
    };
    function searchFilterSortPipelines(pipelines) {
        pipelines = pipelines.filter(filterPipelines).sort(sortPipelines).filter(searchPipelines);
        Filters.set(
            $Filters.map((filter) => {
                if (filter.name === 'Released') {
                    return {
                        name: filter.name,
                        count: pipelines.filter((p) => p.releases.length > 1 && !p.archived).length,
                    };
                }
                if (filter.name === 'Under development') {
                    return {
                        name: filter.name,
                        count: pipelines.filter((p) => p.releases.length === 1 && !p.archived).length,
                    };
                }
                if (filter.name === 'Archived') {
                    return { name: filter.name, count: pipelines.filter((p) => p.archived).length };
                }
                return filter;
            })
        );
        return pipelines;
    }
    SortBy.subscribe(() => {
        filteredPipelines = searchFilterSortPipelines(pipelines);
    });
    CurrentFilter.subscribe(() => {
        filteredPipelines = searchFilterSortPipelines(pipelines);
    });
    SearchQuery.subscribe(() => {
        filteredPipelines = searchFilterSortPipelines(pipelines);
    });

    $: filteredPipelines = searchFilterSortPipelines(pipelines);
    // update counts in CurrentFilter
</script>

<div class="listing d-flex flex-wrap w-100 justify-content-center">
    {#if $DisplayStyle === 'grid'}
        {#if filteredPipelines.length === 0}
            <div class="alert alert-warning" role="alert">
                No pipelines found. Try changing your search query or filters.
            </div>
        {:else}
            {#each filteredPipelines as pipeline (pipeline.name)}
                <PipelineCard {pipeline} />
            {/each}
        {/if}
    {:else if $DisplayStyle === 'table'}
        <table class="table">
            <thead>
                <tr>
                    <th
                        class="text-nowrap sortable"
                        scope="col"
                        data-bs-toggle="tooltip"
                        data-bs-delay="500"
                        title="Sort by name"
                        on:click={() => handleSort('Alphabetical')}
                        ><i
                            class="fa-arrow-up-arrow-down me-2 fa-swap-opacity"
                            class:fa-duotone={$SortBy.startsWith('Alphabetical')}
                            class:fa-regular={!$SortBy.startsWith('Alphabetical')}
                            class:text-muted={!$SortBy.startsWith('Alphabetical')}
                        /> Name</th
                    >
                    <th scope="col">Description</th>
                    <th scope="col">Released</th>
                    <th
                        class="text-end text-nowrap sortable"
                        scope="col"
                        data-bs-toggle="tooltip"
                        data-bs-delay="500"
                        title="Sort by number of stars"
                        on:click={() => handleSort('Stars')}
                        ><i
                            class="fa-arrow-up-arrow-down me-2 fa-swap-opacity"
                            class:fa-duotone={$SortBy.startsWith('Stars')}
                            class:fa-regular={!$SortBy.startsWith('Stars')}
                            class:text-muted={!$SortBy.startsWith('Stars')}
                        /> Stars</th
                    >
                    <th
                        class="text-end text-nowrap sortable"
                        scope="col"
                        data-bs-toggle="tooltip"
                        data-bs-delay="500"
                        title="Sort by date of last release"
                        on:click={() => handleSort('Last release')}
                    >
                        <i
                            class="fa-arrow-up-arrow-down me-2 fa-swap-opacity"
                            class:fa-duotone={$SortBy.startsWith('Last release')}
                            class:fa-regular={!$SortBy.startsWith('Last release')}
                            class:text-muted={!$SortBy.startsWith('Last release')}
                        />Last Release</th
                    >
                </tr>
            </thead>
            <tbody>
                {#each filteredPipelines as pipeline}
                    <tr>
                        <td>
                            <a href={pipeline.html_url} target="_blank" rel="noreferrer">{pipeline.name}</a>
                        </td>
                        <td class="text-small">
                            {pipeline.description}
                        </td>
                        <td class="text-center">
                            {#if pipeline.archived}
                                <i class="fa-solid fa-archive text-info" title="archived" data-bs-toggle="tooltip" />
                            {:else if pipeline.releases.length === 1}
                                <i
                                    class="fa-solid fa-xs fa-wrench text-warning"
                                    title="under development"
                                    data-bs-toggle="tooltip"
                                />
                            {:else if pipeline.releases.length > 1}
                                <i class="fa-solid fa-check text-success" title="released" data-bs-toggle="tooltip" />
                            {/if}
                        </td>
                        <td class="text-end">
                            {pipeline.stargazers_count}
                        </td>
                        <td class="text-end">
                            <a
                                class=""
                                href={'/' +
                                    pipeline.name +
                                    '/' +
                                    (pipeline.releases.length > 1 ? pipeline.releases[0].tag_name : 'dev') +
                                    '/'}
                            >
                                {pipeline.releases.length > 1 ? pipeline.releases[0].tag_name : '-'}
                            </a>
                        </td>
                    </tr>
                {/each}
            </tbody>
        </table>
    {/if}
</div>

<style lang="scss">
    @import '@styles/_variables.scss';
    .table .sortable {
        cursor: pointer;
        &:hover {
            background-color: $secondary-bg-subtle;
        }
        :global([data-bs-theme='dark']) &:hover {
            color: $white;
            background-color: $secondary-bg-subtle-dark;
        }
    }
</style>
