<script lang="ts">
    import ListingTableHeader from '@components/ListingTableHeader.svelte';
    import PipelineCard from '@components/pipeline/PipelineCard.svelte';
    import { CurrentFilter, Filters, SortBy, DisplayStyle, SearchQuery } from '@components/store';
    import { onMount } from 'svelte';

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

    export let filters: { name: string }[] = [{ name: '' }];

    let sortInverse = false;

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
        if ($CurrentFilter.find((f) => f.name === 'Released') && pipeline.releases.length > 1 && !pipeline.archived) {
            return true;
        }
        if (
            $CurrentFilter.find((f) => f.name === 'Under development') &&
            pipeline.releases.length === 1 &&
            !pipeline.archived
        ) {
            return true;
        }
        if ($CurrentFilter.find((f) => f.name === 'Archived') && pipeline.archived === true) {
            return true;
        }
        return false;
    };

    const sortPipelines = (a, b) => {
        sortInverse = $SortBy.endsWith(';inverse');
        if ($SortBy.startsWith('Name')) {
            if (sortInverse) {
                return b.name.localeCompare(a.name);
            } else {
                return a.name.localeCompare(b.name);
            }
        } else if ($SortBy.startsWith('Stars')) {
            if (sortInverse) {
                return a.stargazers_count - b.stargazers_count;
            } else {
                return b.stargazers_count - a.stargazers_count;
            }
        } else if ($SortBy.startsWith('Last release')) {
            // handle case where a pipeline has no releases
            if (a.releases.length === 1 && b.releases.length === 1) {
                if (sortInverse) {
                    return new Date(a.releases[0].published_at) - new Date(b.releases[0].published_at);
                } else {
                    return new Date(b.releases[0].published_at) - new Date(a.releases[0].published_at);
                }
            }
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
            }),
        );
        return pipelines;
    }
    $: filteredPipelines = searchFilterSortPipelines(pipelines);

    onMount(() => {
        console.log(filters);
        CurrentFilter.set(filters);
        SortBy.subscribe(() => {
            filteredPipelines = searchFilterSortPipelines(pipelines);
        });
        CurrentFilter.subscribe(() => {
            filteredPipelines = searchFilterSortPipelines(pipelines);
        });
        SearchQuery.subscribe(() => {
            filteredPipelines = searchFilterSortPipelines(pipelines);
        });
    });
</script>

<div class="listing grid px-2 py-4">
    {#if $DisplayStyle === 'grid'}
        {#if filteredPipelines.length === 0 && $SearchQuery !== ''}
            <div class="g-col-12 g-col-md-8 g-start-md-3">
                <div class="alert alert-secondary text-center" role="alert">
                    No pipelines found. Try changing your search query or filters.
                </div>
            </div>
        {:else}
            {#each filteredPipelines as pipeline (pipeline.name)}
                <div class="g-col-12 g-col-md-6 g-col-xl-4 g-col-xxl-3">
                    <PipelineCard {pipeline} />
                </div>
            {/each}
        {/if}
    {:else}
        <table class="table table-hove table-responsive mx-3">
            <thead>
                <tr>
                    <ListingTableHeader name="Name" />
                    <th scope="col">Description</th>
                    <th scope="col">Released</th>
                    <ListingTableHeader name="Stars" textEnd={true} />
                    <ListingTableHeader name="Last release" title="Sort by date of last release" textEnd={true} />
                </tr>
            </thead>
            <tbody>
                {#each filteredPipelines as pipeline}
                    <tr>
                        <td class=" name p-0">
                            <div class="position-relative p-3">
                                <a
                                    class="stretched-link"
                                    href={'/' + pipeline.name + '/' + pipeline.releases[0].tag_name + '/'}
                                    >{pipeline.name}</a
                                >
                            </div>
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
                            <span>
                                {pipeline.releases.length > 1 ? pipeline.releases[0].tag_name : '-'}
                            </span>
                        </td>
                    </tr>
                {/each}
            </tbody>
        </table>
    {/if}
</div>

<style lang="scss">
</style>
