<script lang="ts">
    import { SearchQuery } from './store.js';

    export let configs: {
        name: string;
        content: string;
        config: object;
    }[] = [];

    const searchconfigs = (config) => {
        if ($SearchQuery === '') {
            return true;
        }
        if (config.meta.name.toLowerCase().includes($SearchQuery.toLowerCase())) {
            return true;
        }
        if (config.content && config.content.toLowerCase().includes($SearchQuery.toLowerCase())) {
            return true;
        }
        return false;
    };

    function searchFilterConfigs(configs) {
        return configs.filter(searchconfigs);
    }
    SearchQuery.subscribe(() => {
        filteredConfigs = searchFilterConfigs(configs);
    });

    $: filteredConfigs = searchFilterConfigs(configs);
    configs.map((config) => {
        config.name = config.name.replace('.md', '');
    });
</script>

<div class="listing d-flex flex-wrap w-100 justify-content-center">
    <!-- {#if $DisplayStyle === 'grid'}
        {#each filteredConfigs as config (config.name)}
            <p>{config.name}</p>
        {/each}
    {:else if $DisplayStyle === 'table'} -->
    <table class="table">
        <thead>
            <tr>
                <th class="name" scope="col">Name</th>
                <th class="description" scope="col">Description</th>
                <th scope="col">Executor</th>
            </tr>
        </thead>
        <tbody>
            {#each filteredConfigs as config}
                <tr>
                    <td class="name">
                        <a href={'/configs/' + config.name}>{@html config.name.replace('_', '_<wbr>')}</a>
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
    <!-- {/if} -->
</div>

<style>
    .name {
        min-width: 15rem;
        word-break: break-word;
    }

    .description {
        min-width: 20rem;
        word-break: break-word;
    }
    .keywords {
        max-width: 35rem;
    }
</style>
