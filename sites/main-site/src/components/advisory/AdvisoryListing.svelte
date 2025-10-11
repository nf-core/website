<script lang="ts">
    import FilterBar from "@components/FilterBar.svelte";
    import AdvisoryCard from "@components/advisory/AdvisoryCard.svelte";
    import { CurrentFilter, SearchQuery } from "@components/store";
    import { onMount } from "svelte";
    import type { CollectionEntry } from "astro:content";
    import { advisoryTypes, advisoryClasses, advisoryIcons } from "./advisoryUtils";
    import { formatAdvisoryType } from "./advisoryUtils";

    interface Props {
        advisories?: CollectionEntry<"advisories">[];
        currentFilters: { name: string }[];
    }

    let { advisories = [], currentFilters }: Props = $props();

    let filteredAdvisories = $state(advisories);

    let sortedAdvisories = $derived(
        [...filteredAdvisories].sort((a, b) => {
            const dateA = a.data.publishedDate?.getTime() ?? 0;
            const dateB = b.data.publishedDate?.getTime() ?? 0;
            return dateB - dateA; // Sort by most recent first
        }),
    );

    const formattedAdvisoryTypes = advisoryTypes.map((type) => ({
        name: type.name,
        displayName: formatAdvisoryType(type.name),
        icon: advisoryIcons[type.name],
        class: advisoryClasses[type.name],
    }));

    const filterByType = (advisories: CollectionEntry<"advisories">) => {
        if ($CurrentFilter.length === 0) return true;
        return $CurrentFilter.some((f) =>
            advisories.data.type.includes(f.name as (typeof advisories.data.type)[number]),
        );
    };

    const searchAdvisories = (advisories: CollectionEntry<"advisories">) => {
        if ($SearchQuery === "") {
            return true;
        }
        // return true if it is in any element of advisories.data
        if (
            Object.values(advisories.data).some((value) => {
                if (typeof value === "string") {
                    return value.toLowerCase().includes($SearchQuery.toLowerCase());
                }
                return false;
            })
        ) {
            return true;
        }
        return false;
    };

    $effect(() => {
        filteredAdvisories = advisories.filter(filterByType).filter(searchAdvisories);
    });

    function hasYearChanged(advisories, idx) {
        if (
            idx === 0 ||
            advisories[idx].data.publishedDate?.getFullYear() !== advisories[idx - 1].data.publishedDate?.getFullYear()
        ) {
            return true;
        }
        return false;
    }

    onMount(() => {
        if (currentFilters.length > 0) {
            CurrentFilter.set(currentFilters);
        }
    });
</script>

<div>
    <FilterBar filter={formattedAdvisoryTypes} displayStyle={[]} sortBy={[]} filterName={() => "Advisory type"}
    ></FilterBar>
    <div class="advisories">
        <div class="d-flex flex-column">
            <div class="mb-3">
                {#each sortedAdvisories as advisories, idx (advisories.id)}
                    {#if hasYearChanged(sortedAdvisories, idx)}
                        <h3 id={"year-" + advisories.data.publishedDate?.getFullYear()}>
                            {advisories.data.publishedDate?.getFullYear()}
                        </h3>
                    {/if}
                    <AdvisoryCard frontmatter={advisories.data} slug={advisories.id} />
                {/each}
            </div>
        </div>
    </div>
</div>

<style lang="scss">
    h3:not(:first-child) {
        margin-top: 0.5rem;
    }
</style>
