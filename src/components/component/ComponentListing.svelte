<script lang="ts">
    import ListingTableHeader from '@components/ListingTableHeader.svelte';
    import ComponentCard from '@components/component/ComponentCard.svelte';
    import { SortBy, DisplayStyle, SearchQuery } from '@components/store';

    export let components: {
        name: string;
        path: string;
        type: string;
        meta: {
            description: string;
            name: string;
        };
        pipelines: {
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
    SortBy.subscribe(() => {
        filteredComponents = searchFilterSortComponents(components);
    });
    SearchQuery.subscribe(() => {
        filteredComponents = searchFilterSortComponents(components);
    });

    $: filteredComponents = searchFilterSortComponents(components);
</script>

<div class={`listing d-flex flex-wrap w-100 justify-content-center ${components[0].type}`}>
    {#if $DisplayStyle === 'grid'}
        {#each filteredComponents as component (component.name)}
            <ComponentCard {component} />
        {/each}
    {:else if $DisplayStyle === 'table'}
        <table class="table">
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
                {#each filteredComponents as component}
                    <tr>
                        <td class="name">
                            <a href={'/' + component.type + 's/' + component.name + '/'}
                                >{@html component.name.replace('_', '_<wbr>')}</a
                            >
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

<style>
    .name {
        min-width: 15rem;
        word-break: break-word;
    }
    .keywords {
        max-width: 35rem;
    }
</style>
