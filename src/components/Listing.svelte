<script lang="ts">
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

    import PipelineCard from './PipelineCard.svelte';
    import { CurrentFilter, SortBy, DisplayStyle, SearchQuery } from './store.js';

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
        if ($CurrentFilter.includes('Released') && pipeline.releases.length > 0 && !pipeline.archived) {
            return true;
        }
        if ($CurrentFilter.includes('Under development') && pipeline.releases.length === 0 && !pipeline.archived) {
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
        return pipelines.filter(filterPipelines).sort(sortPipelines).filter(searchPipelines);
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
</script>

<div class="listing d-flex flex-wrap w-100 justify-content-center">
    {#if $DisplayStyle === 'grid'}
        {#each filteredPipelines as pipeline (pipeline.name)}
            <PipelineCard {pipeline} />
        {/each}
    {:else if $DisplayStyle === 'table'}
        <table class="table">
            <thead>
                <tr>
                    <th scope="col">Name</th>
                    <th scope="col">Description</th>
                    <th scope="col">Status</th>
                    <th scope="col">Stars</th>
                    <th scope="col">Last Release</th>
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
                                : pipeline.releases.length > 0
                                ? 'Released'
                                : 'Under Development'}
                        </td>
                        <td class="text-end">
                            {pipeline.stargazers_count}
                        </td>
                        <td class="text-end">
                            {pipeline.releases.length > 0
                                ? pipeline.releases[pipeline.releases.length - 1].tag_name
                                : '-'}
                        </td>
                    </tr>
                {/each}
            </tbody>
        </table>
    {/if}
</div>

<style>
</style>
