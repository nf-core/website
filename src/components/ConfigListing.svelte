<script lang="ts">
    import ListingTableHeader from '@components/ListingTableHeader.svelte';
    import { SearchQuery, SortBy } from '@components/store';

    export let configs: {
        name: string;
        content: string;
        config: {};
    }[] = [];

    SortBy.set('Name');

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
            if (a.config.executor && b.config.executor) {
                return a.config.executor.localeCompare(b.config.executor) * (invertSort ? -1 : 1);
            } else if (a.config.executor) {
                return -1;
            } else if (b.config.executor) {
                return 1;
            } else {
                return 0;
            }
        }
    };
    function searchFilterSortConfigs(configs) {
        return configs.filter(searchconfigs).sort(sortConfigs);
    }

    SearchQuery.subscribe(() => {
        filteredConfigs = searchFilterSortConfigs(configs);
    });

    SortBy.subscribe(() => {
        filteredConfigs = searchFilterSortConfigs(configs);
    });

    $: filteredConfigs = searchFilterSortConfigs(configs);
    configs.map((config) => {
        config.name = config.name.replace('.md', '');
    });

    const executor_icons = {
        slurm: 'slurm',
        awsbatch: 'aws',
        azurebatch: 'azure',
    };
</script>

<div class="listing d-flex flex-wrap w-100 justify-content-center">
    <table class="table">
        <thead class="text-bg-secondary">
            <tr>
                <ListingTableHeader name="Name" />
                <th class="description" scope="col">Description</th>
                <ListingTableHeader name="Executor" textCenter={true} />
            </tr>
        </thead>
        <tbody>
            {#each filteredConfigs as config}
                <tr class="position-relative">
                    <td class="name">
                        <a class="stretched-link" href={'/configs/' + config.name + '/'}
                            >{@html config.name.replace('_', '_<wbr>')}</a
                        >
                    </td>
                    <td class="description text-small">
                        {config.config.config_profile_description}
                    </td>
                    <td class="text-center">
                        {#if config.config.executor}
                            {#each config.config.executor.split(',') as executor}
                                <span class="badge border border-secondary-subtle text-body fw-normal ms-2">
                                    {executor}
                                </span>
                            {/each}
                        {/if}
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
        vertical-align: middle;
    }
</style>
