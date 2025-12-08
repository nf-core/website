<script lang="ts">
    import FilterBar from "@components/FilterBar.svelte";
    import EventCard from "@components/event/EventCard.svelte";
    import { CurrentFilter, SearchQuery } from "@components/store";
    import { onMount } from "svelte";
    import type { CollectionEntry } from "astro:content";

    interface Props {
        events?: CollectionEntry<"events">[];
        currentFilters: { name: string }[];
        currentEvents: CollectionEntry<"events">[];
    }
    let { events = [], currentFilters, currentEvents = $bindable() }: Props = $props();

    const filterByType = (event: CollectionEntry<"events">) => {
        if ($CurrentFilter.find((f) => f.name === event.data.type)) {
            return true;
        }
        return false;
    };

    const searchEvents = (event: CollectionEntry<"events">) => {
        if ($SearchQuery === "") {
            return true;
        }
        // return true if it is in any element of event.data
        if (
            Object.values(event.data).some((value) => {
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

    let filteredEvents = $derived(events.filter(filterByType).filter(searchEvents));

    function hasRequiredDates(
        event: CollectionEntry<"events">,
    ): event is CollectionEntry<"events"> & { data: { start: Date; end: Date } } {
        return event.data.start !== undefined && event.data.end !== undefined;
    }

    let futureEvents = $derived(
        filteredEvents
            .filter(hasRequiredDates)
            .filter((event) => {
                const today = new Date();
                return event.data.start > today;
            })
            .sort((a, b) => {
                if (a.data.start < b.data.start) {
                    return -1;
                }
                return 1;
            }),
    );

    let pastEvents = $derived(
        filteredEvents
            .filter((event) => {
                const today = new Date();
                return event.data.end && event.data.end < today;
            })
            .sort((a, b) => {
                if (a.data.end && b.data.end && a.data.end < b.data.end) {
                    return 1;
                }
                return -1;
            }),
    );

    let currentEventsFiltered = $derived(
        filteredEvents.filter((event) => {
            const today = new Date();
            return event.data.start && event.data.start < today && event.data.end && event.data.end > today;
        }),
    );

    $effect(() => {
        currentEvents = currentEventsFiltered;
    });

    const event_type_classes = {
        bytesize: "success",
        hackathon: "primary",
        talk: "info",
        training: "warning",
    };
    const event_type_icons = {
        bytesize: "fa-solid fa-apple-core",
        hackathon: "fa-solid fa-laptop-code",
        talk: "fa-solid fa-presentation",
        training: "fa-solid fa-chalkboard-teacher",
    };
    const event_types = Object.keys(event_type_classes).map((type) => {
        return {
            name: type,
            class: event_type_classes[type],
            icon: event_type_icons[type],
        };
    });

    function hasYearChanged(events, idx) {
        if (idx === 0 || events[idx].data.start.getFullYear() !== events[idx - 1].data.start.getFullYear()) {
            return true;
        }
        return false;
    }

    onMount(async () => {
        if (currentFilters.length > 0) {
            CurrentFilter.set(currentFilters);
        }
    });
</script>

<div>
    <FilterBar filter={event_types} displayStyle={[]} sortBy={[]} filterName={() => "Event type"}></FilterBar>
    <div class="events m-auto">
        {#if currentEvents.length > 0}
            <div class="mb-3 col-12">
                <h2><i class="fa-duotone fa-calendar-exclamation me-3"></i>Currently ongoing</h2>
                {#each currentEvents as event (event.id)}
                    <EventCard
                        frontmatter={event.data}
                        slug={event.id}
                        type={event.data.type}
                        time_category="current"
                    />
                {/each}
            </div>
        {/if}
        <div class="mt-5">
            <div class="d-flex flex-column">
                <div class="mb-3">
                    <h2><i class="fa-duotone fa-calendar-day me-3"></i>Upcoming events</h2>
                    {#if futureEvents && futureEvents.length > 0}
                        {#each futureEvents as event (event.id)}
                            <EventCard
                                frontmatter={event.data}
                                slug={event.id}
                                type={event.data.type}
                                time_category="future"
                            />
                        {/each}
                    {:else if !$SearchQuery && $CurrentFilter.length !== 0}
                        <p>Nothing in the calendar at the moment.</p>
                    {/if}
                </div>
                <div class="mb-3">
                    <h2><i class="fa-duotone fa-calendar-check me-3"></i>Past events</h2>
                    {#each pastEvents as event, idx (event.id)}
                        {#if hasYearChanged(pastEvents, idx)}
                            <h3 id={"year-" + event.data.start?.getFullYear()}>{event.data.start?.getFullYear()}</h3>
                        {/if}
                        <EventCard
                            frontmatter={event.data}
                            slug={event.id}
                            type={event.data.type}
                            time_category="past"
                        />
                    {/each}
                </div>
            </div>
        </div>
    </div>
</div>

<style lang="scss">
    .events {
        max-width: 50rem;
    }
</style>
