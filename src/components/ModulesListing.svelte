<script lang="ts">
    import { SortBy, DisplayStyle, SearchQuery } from './store.js';
    import ModuleCard from '@components/ModuleCard.svelte';

    export let modules: {
        name: string;
        meta: {
            description: string;
            name: string;
        };
        pipelines: {
            name: string;
            version: string;
        }[];
    }[] = [];

    const searchModules = (module) => {
        if ($SearchQuery === '') {
            return true;
        }
        if (module.meta.name.toLowerCase().includes($SearchQuery.toLowerCase())) {
            return true;
        }
        if (module.meta.description && module.meta.description.toLowerCase().includes($SearchQuery.toLowerCase())) {
            return true;
        }
        if (
            module.meta.keywords &&
            module.meta.keywords.some((keyword) => keyword.toLowerCase().includes($SearchQuery.toLowerCase()))
        ) {
            return true;
        }
        return false;
    };

    const sortModules = (a, b) => {
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
    function searchFilterSortModules(modules) {
        return modules.sort(sortModules).filter(searchModules);
    }
    SortBy.subscribe(() => {
        filteredModules = searchFilterSortModules(modules);
    });
    SearchQuery.subscribe(() => {
        filteredModules = searchFilterSortModules(modules);
    });

    $: filteredModules = searchFilterSortModules(modules);
</script>

<div class="listing d-flex flex-wrap w-100 justify-content-center">
    {#if $DisplayStyle === 'grid'}
        {#each filteredModules as module (module.name)}
            <ModuleCard {module} />
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
                {#each filteredModules as module}
                    <tr>
                        <td class="name">
                            <a href={'/modules/' + module.name}>{@html module.name.replace('_', '_<wbr>')}</a>
                        </td>
                        <td class="keywords">
                            {#if module.meta.keywords}
                                {#each module.meta.keywords as keyword}
                                    <span class="badge bg-secondary me-1">{keyword}</span>
                                {/each}
                            {:else}
                                -
                            {/if}
                        </td>
                        <td>
                            {module.meta.description}
                        </td>
                        <td class="text-end">
                            {#if module.pipelines}
                                {module.pipelines.length}
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
