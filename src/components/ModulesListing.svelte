<script lang="ts">
    export let modules: {
        name: string;
        meta: {
            description: string;
            name: string;
        };
    }[] = [];
    export let pipelines;
    import { CurrentFilter, SortBy, DisplayStyle, SearchQuery } from './store.js';

    // const searchPipelines = (pipeline) => {
    //     if ($SearchQuery === '') {
    //         return true;
    //     }
    //     if (pipeline.name.toLowerCase().includes($SearchQuery.toLowerCase())) {
    //         return true;
    //     }
    //     if (pipeline.description && pipeline.description.toLowerCase().includes($SearchQuery.toLowerCase())) {
    //         return true;
    //     }
    //     if (pipeline.topics.some((topic) => topic.toLowerCase().includes($SearchQuery.toLowerCase()))) {
    //         return true;
    //     }
    //     return false;
    // };

    // const filterPipelines = (pipeline) => {
    //     if ($CurrentFilter.includes('Released') && pipeline.releases.length > 0 && !pipeline.archived) {
    //         return true;
    //     }
    //     if ($CurrentFilter.includes('Under development') && pipeline.releases.length === 0 && !pipeline.archived) {
    //         return true;
    //     }
    //     if ($CurrentFilter.includes('Archived') && pipeline.archived === true) {
    //         return true;
    //     }
    //     return false;
    // };

    // const sortPipelines = (a, b) => {
    //     if ($SortBy === 'Alphabetical') {
    //         return a.name.localeCompare(b.name);
    //     } else if ($SortBy === 'Stars') {
    //         return b.stargazers_count - a.stargazers_count;
    //     } else if ($SortBy === 'Last release') {
    //         // handle case where a pipeline has no releases
    //         if (a.releases.length === 1) {
    //             return 1;
    //         }
    //         if (b.releases.length === 1) {
    //             return -1;
    //         }
    //         return new Date(b.releases[0].published_at) - new Date(a.releases[0].published_at);
    //     }
    // };
    // function searchFilterSortPipelines(pipelines) {
    //     return pipelines.filter(filterPipelines).sort(sortPipelines).filter(searchPipelines);
    // }
    // SortBy.subscribe(() => {
    //     filteredPipelines = searchFilterSortPipelines(pipelines);
    // });
    // CurrentFilter.subscribe(() => {
    //     filteredPipelines = searchFilterSortPipelines(pipelines);
    // });
    // SearchQuery.subscribe(() => {
    //     filteredPipelines = searchFilterSortPipelines(pipelines);
    // });

    // $: filteredPipelines = searchFilterSortPipelines(pipelines);
    $: filteredModules = modules;

    // count number of pipelines with module in releases
    pipelines.remote_workflows.forEach((pipeline) => {
        const release = pipeline.releases[0];
        if (!release.modules) return;
        release.modules.forEach((module) => {
            const index = modules.findIndex((m) => m.name === module);
            if (index > -1) {
                modules[index].count = modules[index].count ? modules[index].count + 1 : 1;
            }
        });
    });
</script>

<div class="listing d-flex flex-wrap w-100 justify-content-center">
    <!-- {#if $DisplayStyle === 'grid'}
        {#each filteredModules as module (module.name)}
            <PipelineCard {pipeline} />
        {/each}
    {:else if $DisplayStyle === 'table'} -->
    <table class="table">
        <thead>
            <tr>
                <th scope="col">Name</th>
                <th scope="col">Description</th>
                <th class="text-end" scope="col">in # pipelines</th>
                <!-- <th class="text-end" scope="col">Stars</th> -->
                <!-- <th class="text-end" scope="col">Last Release</th> -->
            </tr>
        </thead>
        <tbody>
            {#each filteredModules as module}
                <tr>
                    <td>
                        <a href={module.name}>{module.name}</a>
                    </td>
                    <td>
                        {module.meta.description}
                    </td>
                    <td class="text-end">
                        {#if module.count}
                            {module.count}
                        {:else}
                            -
                        {/if}
                    </td>
                    <!-- <td>
                            {pipeline.archived
                                ? 'Archived'
                                : pipeline.releases.length > 1
                                ? 'Released'
                                : 'Under Development'}
                        </td> -->
                    <!-- <td class="text-end">
                            {pipeline.stargazers_count}
                        </td> -->
                    <!-- <td class="text-end">
                            {pipeline.releases.length > 1
                                ? pipeline.releases[0].tag_name
                                : '-'}
                        </td> -->
                </tr>
            {/each}
        </tbody>
    </table>
    <!-- {/if} -->
</div>

<style>
</style>
