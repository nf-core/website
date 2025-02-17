<script lang="ts">
    import ListingTableHeader from "@components/ListingTableHeader.svelte";
    import { SearchQuery, SortBy } from "@components/store";
    import { type CollectionEntry } from "astro:content";

    interface Props {
        configs?: CollectionEntry<"configs">[];
    }

    let { configs }: Props = $props();
    SortBy.set("Name");

    const searchconfigs = (config) => {
        if ($SearchQuery === "") {
            return true;
        }
        if (config.id.toLowerCase().includes($SearchQuery.toLowerCase())) {
            return true;
        }
        if (config.rendered.metadata.config_profile_description?.toLowerCase().includes($SearchQuery.toLowerCase())) {
            return true;
        }
        if (config.rendered.metadata.executor?.toLowerCase().includes($SearchQuery.toLowerCase())) {
            return true;
        }
        return false;
    };
    let invertSort = false;
    const sortConfigs = (a, b) => {
        invertSort = $SortBy.endsWith(";inverse");
        console.log(invertSort);
        if ($SortBy.startsWith("Name")) {
            return a.id.localeCompare(b.id) * (invertSort ? -1 : 1);
        } else if ($SortBy.startsWith("Executor")) {
            if (a.rendered.metadata.executor && b.rendered.metadata.executor) {
                console.log(a.rendered.metadata.executor, b.rendered.metadata.executor);
                return a.rendered.metadata.executor.localeCompare(b.rendered.metadata.executor) * (invertSort ? -1 : 1);
            } else if (a.rendered.metadata.executor) {
                return -1;
            } else if (b.rendered.metadata.executor) {
                return 1;
            } else {
                return 0;
            }
        }
    };
    function searchFilterSortConfigs(configs) {
        return configs.filter(searchconfigs).sort(sortConfigs);
    }

    let filteredConfigs = $derived(searchFilterSortConfigs(configs));
    $inspect(filteredConfigs[2]);
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
                <tr>
                    <td class="name p-0">
                        <div class="position-relative p-3">
                            <a class="stretched-link" href={"/configs/" + config.id + "/"}
                                >{@html config.id.replaceAll("_", "_<wbr>").replace("conf/", "")}</a
                            >
                        </div>
                    </td>
                    <td class="description text-small">
                        {config.rendered.metadata.config_profile_description}
                    </td>
                    <td class="text-center">
                        {#if config.rendered.metadata.executor}
                            {#each config.rendered.metadata.executor.split(",") as executor}
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
