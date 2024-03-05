<script lang="ts">
    import ListingTableHeader from '@components/ListingTableHeader.svelte';
    import ListingCard from './ListingCard.svelte';
    import { CurrentFilter, Filters, SortBy, DisplayStyle, SearchQuery } from '@components/store';
    import { onMount } from 'svelte';
    import GitHubProfilePictureExtended from '@components/GitHubProfilePictureExtended.svelte';
    import GitHubProfilePicture from '@components/GitHubProfilePicture.svelte';

    import type { CollectionEntry } from 'astro:content';

    export let groups: CollectionEntry<'special-interest-groups'>[] = [];

    let filteredGroups = groups;
    console.log(filteredGroups);
</script>

<div class="listing px-2 py-4">
    {#if $DisplayStyle === 'grid'}
        <div class="grid">
            {#if filteredGroups.length === 0 && $SearchQuery !== ''}
                <div class="g-col-12 g-col-md-8 g-start-md-3">
                    <div class="alert alert-secondary text-center" role="alert">
                        No Groups found. Try changing your search query or filters.
                    </div>
                </div>
            {:else}
                {#each filteredGroups as group (group.id)}
                    <div class="g-col-12 g-col-lg-6 g-col-xl-6 g-col-xxl-4 g-col-xxxxl-2">
                        <ListingCard footer={true}>
                            <a slot="card-header" href={'/special-interest-groups/' + group.slug} class="success"
                                >{group.data.title}</a
                            >
                            <p slot="card-body">{group.data.subtitle}</p>
                            <div slot="card-footer" class="d-flex align-items-start">
                                <div class="pipeline-badges small flex-grow-1 border-end me-3 h-100">
                                    {#if group.data.pipelines}
                                        <p class="text-muted small mb-1">Pipelines:</p>
                                        {#each group.data.pipelines as pipeline}
                                            <span class={`badge me-2 pipeline-badge small`}>{pipeline}</span>
                                        {/each}
                                    {/if}
                                </div>
                                <div class="leads w-75 small text-end">
                                    {#if group.data.leads}
                                        <p class="text-muted small mb-2">Group leads:</p>
                                        {#each group.data.leads as lead}
                                            {#if typeof lead === 'string'}
                                                <GitHubProfilePictureExtended username={lead} size={25} />
                                            {:else}
                                                <GitHubProfilePictureExtended
                                                    username={Object.keys(lead)[0]}
                                                    size={25}
                                                    wrapperClasses=" justify-content-end"
                                                >
                                                    {Object.values(lead)[0]}
                                                </GitHubProfilePictureExtended>
                                            {/if}
                                        {/each}
                                    {/if}
                                </div>
                            </div></ListingCard
                        >
                    </div>
                {/each}
            {/if}
        </div>
    {:else}
        <table class="table table-hove table-responsive">
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
                                <a class="stretched-link" href={'/special-interest-groups' + group.slug + '/'}
                                    >{group.data.title}</a
                                >
                            </div>
                        </td>
                        <td class="text-small">
                            {group.data.subtitle}
                        </td>
                        <td class="pipeline-badges small">
                            {#each group.data.pipelines ?? [] as pipeline}
                                <span class={`badge me-2 pipeline-badge`}>{pipeline}</span>
                            {/each}
                        </td>
                        <td>
                            {#each group.data.leads ?? [] as lead}
                                {#if typeof lead === 'string'}
                                    <GitHubProfilePictureExtended username={lead} size={25} />
                                {:else}
                                    <GitHubProfilePictureExtended username={Object.keys(lead)[0]} size={25}>
                                        {Object.values(lead)[0]}
                                    </GitHubProfilePictureExtended>
                                {/if}
                            {/each}
                        </td>
                    </tr>
                {/each}
            </tbody>
        </table>
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
