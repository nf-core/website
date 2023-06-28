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
        if ($SortBy === 'Alphabetical') {
            return a.name.localeCompare(b.name);
        } else if ($SortBy === 'Stars') {
            return b.stargazers_count - a.stargazers_count;
        } else if ($SortBy === 'Last release') {
            // handle case where a pipeline has no releases
            if (a.releases.length === 1) {
                return 1;
            }
            if (b.releases.length === 1) {
                return -1;
            }
            return new Date(b.releases[0].published_at) - new Date(a.releases[0].published_at);
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
        console.log($Filters);
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
                    <th scope="col">Name</th>
                    <th scope="col">Description</th>
                    <th scope="col">Status</th>
                    <th class="text-end" scope="col">Stars</th>
                    <th class="text-end" scope="col">Last Release</th>
                </tr>
            </thead>
            <tbody>
                {#each filteredPipelines as pipeline}
                    <tr>
                        <td>
                            <a href={pipeline.html_url} target="_blank" rel="noreferrer">{pipeline.name}</a>
                        </td>
                        <td>
                            {pipeline.description}
                        </td>
                        <td>
                            {pipeline.archived
                                ? 'Archived'
                                : pipeline.releases.length > 1
                                ? 'Released'
                                : 'Under Development'}
                        </td>
                        <td class="text-end">
                            {pipeline.stargazers_count}
                        </td>
                        <td class="text-end">
                            {pipeline.releases.length > 1 ? pipeline.releases[0].tag_name : '-'}
                        </td>
                    </tr>
                {/each}
            </tbody>
        </table>
    {/if}
</div>

<style>
</style>
