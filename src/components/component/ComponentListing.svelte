<script lang="ts">
    import { SortBy, DisplayStyle, SearchQuery } from '@components/store';
    import ComponentCard from '@components/component/ComponentCard.svelte';

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

    const sortComponents = (a, b) => {
        if ($SortBy === 'Alphabetical') {
            return a.name.localeCompare(b.name);
        } else if ($SortBy === '# Pipeline integrations') {
            if (a.pipelines && b.pipelines) {
                return b.pipelines.length - a.pipelines.length;
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

<div class="listing d-flex flex-wrap w-100 justify-content-center">
    {#if $DisplayStyle === 'grid'}
        {#each filteredComponents as component (component.name)}
            <ComponentCard {component} />
        {/each}
    {:else if $DisplayStyle === 'table'}
        <table class="table">
            <thead>
                <tr>
                    <th class="name" scope="col">Name</th>
                    <th class="keywords" scope="col">Keywords</th>
                    <th scope="col">Description</th>
                    <th class="text-end" scope="col">in # pipelines</th>
                </tr>
            </thead>
            <tbody>
                {#each filteredComponents as component}
                    <tr>
                        <td class="name">
                            <a href={'/' + component.type + 's/' + component.name}
                                >{@html component.name.replace('_', '_<wbr>')}</a
                            >
                        </td>
                        <td class="keywords">
                            {#if component.meta.keywords}
                                {#each component.meta.keywords as keyword}
                                    <span class="badge bg-secondary me-1">{keyword}</span>
                                {/each}
                            {:else}
                                -
                            {/if}
                        </td>
                        <td>
                            {component.meta.description}
                        </td>
                        <td class="text-end">
                            {#if component.pipelines}
                                {component.pipelines.length}
                            {:else}
                                -
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
