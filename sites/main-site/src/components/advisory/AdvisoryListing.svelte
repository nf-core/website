<script lang="ts">
    import FilterBar from "@components/FilterBar.svelte";
    import AdvisoryCard from "@components/advisory/AdvisoryCard.svelte";
    import { CurrentFilter, SearchQuery } from "@components/store";
    import { onMount } from "svelte";
    import type { CollectionEntry } from "astro:content";

    interface Props {
        advisories?: CollectionEntry<"advisories">[];
        currentFilters: { name: string }[];
        currentAdvisory: CollectionEntry<"advisories">[];
    }

    let { advisories = [], currentFilters, currentAdvisory = $bindable() }: Props = $props();

    let filteredAdvisorys = $state(advisories);
    const filterByType = (advisories: CollectionEntry<"advisories">) => {
        if ($CurrentFilter.find((f) => f.name === advisories.data.type)) {
            return true;
        }
        return false;
    };

    const searchAdvisorys = (advisories: CollectionEntry<"advisories">) => {
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

    function hasRequiredDates(
        advisories: CollectionEntry<"advisories">,
    ): advisories is CollectionEntry<"advisories"> & { data: { start: Date; end: Date } } {
        return advisories.data.start !== undefined && advisories.data.end !== undefined;
    }

    let futureAdvisorys = $derived(
        filteredAdvisorys
            .filter(hasRequiredDates)
            .filter((advisories) => {
                const today = new Date();
                return advisories.data.start > today;
            })
            .sort((a, b) => {
                if (a.data.start < b.data.start) {
                    return -1;
                }
                return 1;
            }),
    );

    let pastAdvisorys = $derived(
        filteredAdvisorys
            .filter((advisories) => {
                const today = new Date();
                return advisories.data.end && advisories.data.end < today;
            })
            .sort((a, b) => {
                if (a.data.end && b.data.end && a.data.end < b.data.end) {
                    return 1;
                }
                return -1;
            }),
    );

    let currentAdvisorysFiltered = $derived(
        filteredAdvisorys.filter((advisories) => {
            const today = new Date();
            return advisories.data.start && advisories.data.start < today && advisories.data.end && advisories.data.end > today;
        }),
    );

    $effect(() => {
        filteredAdvisorys = advisories.filter(filterByType).filter(searchAdvisorys);
    });

    $effect(() => {
        currentAdvisorys = currentAdvisorysFiltered;
    });

    const advisories_type_classes = {
        bytesize: "success",
        hackathon: "primary",
        talk: "info",
        training: "warning",
    };
    const advisories_type_icons = {
        bytesize: "fa-solid fa-apple-core",
        hackathon: "fa-solid fa-laptop-code",
        talk: "fa-solid fa-presentation",
        training: "fa-solid fa-chalkboard-teacher",
    };
    const advisories_types = Object.keys(advisories_type_classes).map((type) => {
        return {
            name: type,
            class: advisories_type_classes[type],
            icon: advisories_type_icons[type],
        };
    });

    function hasYearChanged(advisories, idx) {
        if (idx === 0 || advisories[idx].data.start.getFullYear() !== advisories[idx - 1].data.start.getFullYear()) {
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
        {#if currentAdvisorys.length > 0}
            <div class="mb-3 col-12">
                <h2><i class="fa-duotone fa-calendar-exclamation me-3"></i>Currently ongoing</h2>
                {#each currentAdvisorys as advisories (advisories.id)}
                    <AdvisoryCard
                        frontmatter={advisories.data}
                        slug={advisories.id}
                        type={advisories.data.type}
                        time_category="current"
                    />
                {/each}
            </div>
        {/if}
        <div class="mt-5">
            <div class="d-flex flex-column">
                <div class="mb-3">
                    <h2><i class="fa-duotone fa-calendar-day me-3"></i>Upcoming advisories</h2>
                    {#if futureAdvisorys && futureAdvisorys.length > 0}
                        {#each futureAdvisorys as advisories (advisories.id)}
                            <AdvisoryCard
                                frontmatter={advisories.data}
                                slug={advisories.id}
                                type={advisories.data.type}
                                time_category="future"
                            />
                        {/each}
                    {:else if !$SearchQuery && $CurrentFilter.length !== 0}
                        <p>Nothing in the calendar at the moment.</p>
                    {/if}
                </div>
                <div class="mb-3">
                    <h2><i class="fa-duotone fa-calendar-check me-3"></i>Past advisories</h2>
                    {#each pastAdvisorys as advisories, idx (advisories.id)}
                        {#if hasYearChanged(pastAdvisorys, idx)}
                            <h3 id={"year-" + advisories.data.start?.getFullYear()}>{advisories.data.start?.getFullYear()}</h3>
                        {/if}
                        <AdvisoryCard
                            frontmatter={advisories.data}
                            slug={advisories.id}
                            type={advisories.data.type}
                            time_category="past"
                        />
                    {/each}
                </div>
            </div>
        </div>
    </div>
</div>
