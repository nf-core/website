<script lang="ts">
    import FilterBar from "@components/FilterBar.svelte";
    import AdvisoryCard from "@components/advisory/AdvisoryCard.svelte";
    import { CurrentFilter, SearchQuery } from "@components/store";
    import { onMount } from "svelte";
    import type { CollectionEntry } from "astro:content";
    import { advisories_types } from "./advisoryTypes";

    interface Props {
        advisories?: CollectionEntry<"advisories">[];
        currentFilters: { name: string }[];
    }

    let { advisories = [], currentFilters }: Props = $props();

    let filteredAdvisories = $state(advisories);
    let currentAdvisories = $derived(
        filteredAdvisories.filter((advisories) => {
            // Calculate date 3 months ago
            const threeMonthsAgo = new Date();
            threeMonthsAgo.setMonth(threeMonthsAgo.getMonth() - 3);
            const threeMonthsAgoTime = threeMonthsAgo.getTime();

            // Check if advisory was published within the last 3 months
            const publishedDateUnix = advisories.data.publishedDate?.getTime();
            return publishedDateUnix && publishedDateUnix >= threeMonthsAgoTime;
        }),
    );

    let pastAdvisories = $derived(
        filteredAdvisories
            .filter((advisories) => {
                // Calculate date 3 months ago
                const threeMonthsAgo = new Date();
                threeMonthsAgo.setMonth(threeMonthsAgo.getMonth() - 3);
                const threeMonthsAgoTime = threeMonthsAgo.getTime();

                // Check if advisory was published more than 3 months ago
                const publishedDateUnix = advisories.data.publishedDate?.getTime();
                return !publishedDateUnix || publishedDateUnix < threeMonthsAgoTime;
            })
            .sort((a, b) => {
                const dateA = a.data.publishedDate?.getTime() ?? 0;
                const dateB = b.data.publishedDate?.getTime() ?? 0;
                return dateB - dateA; // Sort by most recent first
            }),
    );

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
    <FilterBar filter={advisories_types} displayStyle={[]} sortBy={[]} filterName={() => "Advisory type"}></FilterBar>
    <div class="advisories">
        {#if currentAdvisories.length > 0}
            <div class="mb-3 col-12">
                <h2><i class="fa-duotone fa-calendar-exclamation me-3"></i>Recent advisories</h2>
                {#each currentAdvisories as advisories (advisories.id)}
                    <AdvisoryCard frontmatter={advisories.data} slug={advisories.id} time_category="current" />
                {/each}
            </div>
        {/if}
        <div class="mt-5">
            <div class="d-flex flex-column">
                <div class="mb-3">
                    <h2><i class="fa-duotone fa-calendar-check me-3"></i>Past advisories</h2>
                    {#each pastAdvisories as advisories, idx (advisories.id)}
                        {#if hasYearChanged(pastAdvisories, idx)}
                            <h3 id={"year-" + advisories.data.publishedDate?.getFullYear()}>
                                {advisories.data.publishedDate?.getFullYear()}
                            </h3>
                        {/if}
                        <AdvisoryCard frontmatter={advisories.data} slug={advisories.id} time_category="past" />
                    {/each}
                </div>
            </div>
        </div>
    </div>
</div>
