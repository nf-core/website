<script lang="ts">
    import { config } from 'dotenv';
    import { DisplayStyle, SearchQuery } from './store.js';
    import ComponentCard from '@components/ComponentCard.svelte';

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
                <th scope="col">Description</th>
                <th scope="col">URL</th>
            </tr>
        </thead>
        <tbody>
            {#each filteredConfigs as config}
                <tr>
                    <td class="name">
                        <a href={config.name}>{@html config.name.replace('_', '_<wbr>')}</a>
                    </td>
                    <td class="description">
                        {config.config?.params?.config_profile_description
                            ? config.config.params.config_profile_description
                            : ''}
                    </td>
                    <td class="url">
                        {config.config?.params?.config_profile_url
                            ? `<a href=${config.config.params.config_profile_url}>${config.config.params.config_profile_url}</a>`
                            : ''}
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
    .keywords {
        max-width: 35rem;
    }
</style>
