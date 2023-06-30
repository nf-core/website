<script lang="ts">
    import { SearchQuery, SortBy } from '@components/store';
    import ListingTableHeader from '@components/ListingTableHeader.svelte';

    export let configs: {
        name: string;
        content: string;
        config: {};
    }[] = [];

    const searchconfigs = (config) => {
        if ($SearchQuery === '') {
            return true;
        }
        if (config.name.toLowerCase().includes($SearchQuery.toLowerCase())) {
            return true;
        }
        if (config.config_profile_description?.toLowerCase().includes($SearchQuery.toLowerCase())) {
            return true;
        }
        if (config.executor?.toLowerCase().includes($SearchQuery.toLowerCase())) {
            return true;
        }
        return false;
    };
    let invertSort = false;
    const sortConfigs = (a, b) => {
        invertSort = $SortBy.endsWith(';inverse');
        if ($SortBy.startsWith('Name')) {
            return a.name.localeCompare(b.name) * (invertSort ? -1 : 1);
        } else if ($SortBy.startsWith('Executor')) {
            return a.config.executor.localeCompare(b.config.executor) * (invertSort ? -1 : 1);
        }
    };
    function searchFilterSortConfigs(configs) {
        return configs.filter(searchconfigs).sort(sortConfigs);
    }

    SearchQuery.subscribe(() => {
        filteredConfigs = searchFilterSortConfigs(configs);
    });

    $: filteredConfigs = searchFilterSortConfigs(configs);
    configs.map((config) => {
        config.name = config.name.replace('.md', '');
    });
</script>

<div class="listing d-flex flex-wrap w-100 justify-content-center">
    <table class="table">
        <thead class="text-bg-secondary">
            <tr>
                <ListingTableHeader name="Name" />
                <th class="description" scope="col">Description</th>
                <ListingTableHeader name="Executor" />
            </tr>
        </thead>
        <tbody>
            {#each filteredConfigs as config}
                <tr>
                    <td class="name">
                        <a href={'/configs/' + config.name + '/'}>{@html config.name.replace('_', '_<wbr>')}</a>
                    </td>
                    <td class="description w-50">
                        {config.config.config_profile_description}
                    </td>
                    <td>
                        {config.config.executor}
                    </td>
                </tr>
            {/each}
        </tbody>
    </table>
</div>

<style lang="scss">
    .name {
        min-width: 15rem;
        word-break: break-word;
    }

    .description {
        min-width: 20rem;
        word-break: break-word;
    }
</style>
