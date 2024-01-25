<script lang="ts">
    import ListingTableHeader from '@components/ListingTableHeader.svelte';
    import PaginationNav from '@components/PaginationNav.svelte';
    import ComponentCard from '@components/component/ComponentCard.svelte';
    import { SortBy, DisplayStyle, SearchQuery, currentPage } from '@components/store';

    export let components: {
        name: string;
        path: string;
        type: string;
        meta: {
            name: string;
            description: string;
            keywords: string[];
            components?: string[];
            input: {}[];
            output: {}[];
        };
        pipelines?: {
            name: string;
            version: string;
        }[];
    }[] = [];

    const searchComponents = (component) => {
        if ($SearchQuery === '') {
            return true;
        }
        if (component.meta.name.toLowerCase().includes($SearchQuery.toLowerCase())) {
            return true;
        }
        if (
            component.meta.description &&
            component.meta.description.toLowerCase().includes($SearchQuery.toLowerCase())
        ) {
            return true;
        }
        if (
            component.meta.keywords &&
            component.meta.keywords.some((keyword) => keyword.toLowerCase().includes($SearchQuery.toLowerCase()))
        ) {
            return true;
        }
        return false;
    };
    let invertSort = false;
    const sortComponents = (a, b) => {
        invertSort = $SortBy.endsWith(';inverse');
        if ($SortBy.startsWith('Name')) {
            return a.name.localeCompare(b.name) * (invertSort ? -1 : 1);
        } else if ($SortBy.startsWith('# Pipeline integrations')) {
            if (a.pipelines && b.pipelines) {
                return (b.pipelines.length - a.pipelines.length) * (invertSort ? -1 : 1);
            } else if (a.pipelines) {
                return -1;
            } else if (b.pipelines) {
                return 1;
            } else {
                return 0;
            }
        }
    };
    function searchFilterSortComponents(components) {
        return components.sort(sortComponents).filter(searchComponents);
    }

    $: filteredComponents = searchFilterSortComponents(components) || [];

    let pageSize: number = 12;
    let lastPage = Math.ceil(components.length / pageSize);
    const updatePageSize = () => {
        pageSize = $DisplayStyle === 'grid' ? 12 : 25;
        let currentComponents = filteredComponents || components;
        lastPage = Math.ceil(currentComponents.length / pageSize);
    };
    updatePageSize();

    $: paginatedItems = filteredComponents.slice(($currentPage - 1) * pageSize, $currentPage * pageSize);

    SortBy.subscribe(() => {
        filteredComponents = searchFilterSortComponents(components);
    });
    SearchQuery.subscribe(() => {
        filteredComponents = searchFilterSortComponents(components);
        updatePageSize();
    });
    DisplayStyle.subscribe(() => {
        updatePageSize();
    });
</script>

<div class={`listing grid px-2 py-4 ${components[0].type}`}>
    {#if $DisplayStyle === 'grid'}
        {#each paginatedItems as component (component.name)}
            <div class="g-col-12 g-col-md-6 g-col-xl-4 g-col-xxl-3">
                <ComponentCard {component} />
            </div>
        {/each}
    {:else}
        <table class="table table-responsive mx-3">
            <thead>
                <tr>
                    <ListingTableHeader name="Name" />
                    <th scope="col">Description</th>
                    <th class="keywords" scope="col">Keywords</th>
                    <ListingTableHeader
                        name="# Pipeline integrations"
                        title={'Sort by number of pipelines with ' + components[0].type}
                        textEnd={true}
                    />
                </tr>
            </thead>
            <tbody>
                {#each paginatedItems as component (component.name)}
                    <tr>
                        <td class=" name p-0">
                            <div class="position-relative p-3">
                                <a class="stretched-link" href={'/' + component.type + 's/' + component.name + '/'}
                                    >{@html component.name.replace('_', '_<wbr>')}</a
                                >
                            </div>
                        </td>
                        <td class="text-small">
                            {component.meta.description}
                        </td>
                        <td class="topics">
                            <!-- {#if component.meta.keywords} -->
                            {#each component.meta.keywords as keyword}
                                <span class={`badge me-2 ${component.type}-topic`}>{keyword}</span>
                            {/each}
                            <!-- {/if} -->
                        </td>
                        <td class="text-end">
                            {#if component.pipelines}
                                {component.pipelines.length}
                            {/if}
                        </td>
                    </tr>
                {/each}
            </tbody>
        </table>
    {/if}
</div>
{#if lastPage > 0}
    <PaginationNav {lastPage} />
{/if}

<style>
    .name {
        min-width: 15rem;
        word-break: break-word;
    }
    .keywords {
        max-width: 35rem;
    }
</style>
