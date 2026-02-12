<script lang="ts">
    import ListingTableHeader from "@components/ListingTableHeader.svelte";
    import ListingCard from "./ListingCard.svelte";
    import { DisplayStyle, SearchQuery } from "@components/store";

    import type { CollectionEntry } from "astro:content";

    interface Props {
        groups?: CollectionEntry<"special-interest-groups">[];
    }

    let { groups = [] }: Props = $props();

    let filteredGroups = groups;

    const getLeadName = (lead: string | object): string => {
        return typeof lead === "string" ? lead : Object.keys(lead)[0];
    };
</script>

<div class="listing px-2 py-4">
    {#if $DisplayStyle === "grid"}
        <div class="grid">
            {#if filteredGroups.length === 0 && $SearchQuery !== ""}
                <div class="g-col-12 g-col-md-8 g-start-md-3">
                    <div class="alert alert-secondary text-center" role="alert">
                        No Groups found. Try changing your search query or filters.
                    </div>
                </div>
            {:else}
                {#each filteredGroups as group (group.id)}
                    <div class="g-col-12 g-col-lg-6 g-col-xl-6 g-col-xxl-4 g-col-xxxxl-2">
                        <ListingCard footer={true}>
                            {#snippet cardHeader()}
                                <a href={"/special-interest-groups/" + group.id} class="success">
                                    {group.data.groupName}
                                </a>
                            {/snippet}

                            {#snippet cardBody()}
                                <p>{group.data.subtitle}</p>
                            {/snippet}

                            {#snippet cardFooter()}
                                <div class="grid align-content-start">
                                    <div class="pipeline-badges small g-col-8">
                                        {#if group.data.pipelines && group.data.pipelines.length > 0}
                                            <p class="text-muted small mb-1">Pipelines:</p>
                                            {#each group.data.pipelines as pipeline}
                                                <span class={`badge me-2 pipeline-badge small`}>{pipeline}</span>
                                            {/each}
                                        {/if}
                                    </div>
                                    <div class="small g-col-4">
                                        {#if group.data.leads}
                                            <p class="text-muted small mb-2">Group leads:</p>
                                            <div class="leads d-flex w-100 h-100 flex-wrap align-content-start gap-1">
                                                {#each group.data.leads as lead}
                                                    <a
                                                        href={`https://github.com/${getLeadName(lead)}`}
                                                        target="_blank"
                                                        rel="noopener noreferrer"
                                                        class="badge text-bg-dark text-decoration-none"
                                                        data-bs-toggle="tooltip"
                                                        title={getLeadName(lead)}
                                                    >
                                                        @{getLeadName(lead)}
                                                    </a>
                                                {/each}
                                            </div>
                                        {/if}
                                    </div>
                                </div>
                            {/snippet}
                        </ListingCard>
                    </div>
                {/each}
            {/if}
        </div>
    {:else}
        <div class="table-responsive">
            <table class="table table-hover">
                <thead>
                    <tr>
                        <ListingTableHeader name="Name" />
                        <th scope="col">Description</th>
                        <th scope="col">Included pipelines</th>
                        <th scope="col">Leads</th>
                    </tr>
                </thead>
                <tbody>
                    {#each filteredGroups as group (group.id)}
                        <tr>
                            <td class=" name p-0">
                                <div class="position-relative p-3">
                                    <a class="stretched-link" href={"/special-interest-groups/" + group.id}
                                        >{group.data.groupName}</a
                                    >
                                </div>
                            </td>
                            <td class="text-small">
                                {group.data.subtitle}
                            </td>
                            <td class="pipeline-badges small col-3">
                                {#each group.data.pipelines ?? [] as pipeline}
                                    <span class={`badge me-2 pipeline-badge`}>{pipeline}</span>
                                {/each}
                            </td>
                            <td class="">
                                <div class="d-flex flex-wrap gap-1">
                                    {#each group.data.leads ?? [] as lead}
                                        <a
                                            href={`https://github.com/${getLeadName(lead)}`}
                                            target="_blank"
                                            rel="noopener noreferrer"
                                            class="badge text-bg-dark text-decoration-none"
                                            data-bs-toggle="tooltip"
                                            title={getLeadName(lead)}
                                        >
                                            @{getLeadName(lead)}
                                        </a>
                                    {/each}
                                </div>
                            </td>
                        </tr>
                    {/each}
                </tbody>
            </table>
        </div>
    {/if}
</div>

<style lang="scss">
    .leads {
        width: fit-content;
    }
    :global(h2) a {
        color: var(--bs-success);
    }
</style>
