<script lang="ts">
    import ListingTableHeader from "@components/ListingTableHeader.svelte";
    import PipelineCard from "@components/pipeline/PipelineCard.svelte";
    import { CurrentFilter, SortBy, DisplayStyle, SearchQuery } from "@components/store";
    import { onMount } from "svelte";

    interface Pipeline {
        name: string;
        description: string;
        stargazers_count: number;
        topics: string[];
        releases: {
            published_at: string;
            tag_name: string;
        }[];
        archived: boolean;
    }

    let { pipelines, filters } = $props();

    onMount(() => {
        CurrentFilter.set(filters);
    });

    let currentFilter = $derived($CurrentFilter || filters);
    let displayStyle = $derived($DisplayStyle || "grid");
    let searchQuery = $derived($SearchQuery || "");
    let sortBy = $derived($SortBy || "Last release");
    let sortInverse = $derived(sortBy.endsWith(";inverse"));

    const searchPipelines = (pipeline: Pipeline) => {
        if (searchQuery === "") {
            return true;
        }
        if (pipeline.name.toLowerCase().includes(searchQuery.toLowerCase())) {
            return true;
        }
        if (pipeline.description && pipeline.description.toLowerCase().includes(searchQuery.toLowerCase())) {
            return true;
        }
        if (pipeline.topics.some((topic) => topic.toLowerCase().includes(searchQuery.toLowerCase()))) {
            return true;
        }
        return false;
    };

    const filterPipelines = (pipeline: Pipeline) => {
        if (currentFilter.find((f) => f.name === "Released") && pipeline.releases.length > 1 && !pipeline.archived) {
            return true;
        }
        if (
            currentFilter.find((f) => f.name === "Under development") &&
            pipeline.releases.length === 1 &&
            !pipeline.archived
        ) {
            return true;
        }
        if (currentFilter.find((f) => f.name === "Archived") && pipeline.archived === true) {
            return true;
        }
        return false;
    };

    const sortPipelines = (a: Pipeline, b: Pipeline) => {
        if (sortBy.startsWith("Name")) {
            if (sortInverse) {
                return b.name.localeCompare(a.name);
            } else {
                return a.name.localeCompare(b.name);
            }
        } else if (sortBy.startsWith("Stars")) {
            if (sortInverse) {
                return a.stargazers_count - b.stargazers_count;
            } else {
                return b.stargazers_count - a.stargazers_count;
            }
        } else if (sortBy.startsWith("Last release")) {
            if (a.releases.length === 1 && b.releases.length === 1) {
                if (sortInverse) {
                    return (
                        new Date(a.releases[0].published_at).getTime() - new Date(b.releases[0].published_at).getTime()
                    );
                } else {
                    return (
                        new Date(b.releases[0].published_at).getTime() - new Date(a.releases[0].published_at).getTime()
                    );
                }
            }
            if (a.releases.length === 1) {
                return 1 * (sortInverse ? -1 : 1);
            }
            if (b.releases.length === 1) {
                return -1 * (sortInverse ? -1 : 1);
            }

            if (sortInverse) {
                return new Date(a.releases[0].published_at).getTime() - new Date(b.releases[0].published_at).getTime();
            } else {
                return new Date(b.releases[0].published_at).getTime() - new Date(a.releases[0].published_at).getTime();
            }
        }
        return 0;
    };
    let filteredPipelines = $derived(pipelines.filter(filterPipelines).filter(searchPipelines).sort(sortPipelines));
</script>

<div class="listing px-2 py-4">
    {#if displayStyle === "grid"}
        <div class="grid">
            {#if filteredPipelines.length === 0 && searchQuery !== ""}
                <div class="g-col-12 g-col-md-8 g-start-md-3">
                    <div class="alert alert-secondary text-center" role="alert">
                        No pipelines found. Try changing your search query or filters.
                    </div>
                </div>
            {:else}
                {#each filteredPipelines as pipeline (pipeline.name)}
                    <div class="g-col-12 g-col-md-6 g-col-xl-6 g-col-xxl-4 g-col-xxxl-3 g-col-xxxxl-2">
                        <PipelineCard {pipeline} />
                    </div>
                {/each}
            {/if}
        </div>
    {:else}
        <table class="table table-hover table-responsive">
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
                        <td class="name p-0">
                            <div class="position-relative p-3">
                                <a
                                    class="stretched-link"
                                    href={"/" + pipeline.name + "/" + pipeline.releases[0].tag_name + "/"}
                                    >{pipeline.name}</a
                                >
                            </div>
                        </td>
                        <td class="text-small">
                            {pipeline.description}
                        </td>
                        <td class="text-center">
                            {#if pipeline.archived}
                                <i class="fa-solid fa-archive text-info" title="archived" data-bs-toggle="tooltip"></i>
                            {:else if pipeline.releases.length === 1}
                                <i
                                    class="fa-solid fa-xs fa-wrench text-warning"
                                    title="under development"
                                    data-bs-toggle="tooltip"
                                ></i>
                            {:else if pipeline.releases.length > 1}
                                <i class="fa-solid fa-check text-success" title="released" data-bs-toggle="tooltip"></i>
                            {/if}
                        </td>
                        <td class="text-end">
                            {pipeline.stargazers_count}
                        </td>
                        <td class="text-end">
                            <span>
                                {pipeline.releases.length > 1 ? pipeline.releases[0].tag_name : "-"}
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
