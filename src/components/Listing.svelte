<script>
    export let pipelines = [];
import PipelineCard from './PipelineCard.svelte';
import { CurrentFilter, SortBy, DisplayStyle } from './filter.js';

import { Card, CardHeader, CardBody, CardFooter} from 'sveltestrap';
let displayStyle = 'grid'

const filterPipelines = (pipeline) => {
        if ($CurrentFilter.includes('released') && pipeline.releases.length > 0 && pipeline.archived !== true) {
            return true;
        }
        if ($CurrentFilter.includes('under_development') && pipeline.releases.length === 0) {
            return true;
        }
        if ($CurrentFilter.includes('archived') && pipeline.archived === true) {
            return true;
        }
        return false;
    };

    const sortPipelines = (a, b) => {
        if ($SortBy === 'alphabetical') {
            return a.name.localeCompare(b.name);
        } else if ($SortBy === 'stars') {
            return b.stargazers_count - a.stargazers_count;
        } else if ($SortBy === 'last_release') {
            // handle case where pipeline has no releases
            if (a.releases.length === 0) {
                return 1;
            }
            if (b.releases.length === 0) {
                return -1;
            }
            return (
                new Date(b.releases[b.releases.length - 1].published_at) -
                new Date(a.releases[a.releases.length - 1].published_at)
            );
        }
    };

    $: filteredPipelines = pipelines.filter(filterPipelines).sort(sortPipelines);
// A react component to display a grid of pipelines
</script>

<div class="listing d-flex flex-wrap w-100 justify-content-center">
    {#if displayStyle=== 'grid'}
        {#each pipelines as pipeline}
            <!-- <Card class="w-25 m-2">
                <CardHeader>
                    {pipeline.name}
                </CardHeader>
                <CardBody>
                   {pipeline.description}
                </CardBody>

            </Card> -->
 <PipelineCard key={pipeline.name} pipeline={pipeline} />
        {/each}
    {:else}
        {#each filteredPipelines as pipeline}
            <!-- <PipelineCard key={pipeline.name} pipeline={pipeline} /> -->
        {/each}
    {/if}
</div>
<style>
    .card {
        width: 300px;
    }
</style>
