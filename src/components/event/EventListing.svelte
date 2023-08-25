<script lang="ts">
    import FilterBar from '@components/FilterBar.svelte';
    import EventCard from '@components/event/EventCard.svelte';
    import { CurrentFilter, SearchQuery, EventIsOngoing } from '@components/store';
    import { onMount } from 'svelte';

    export let events = [];
    export let currentFilters: { name: string }[];
    export let currentEvents;

    let filteredEvents = events;
    const filterByType = (event) => {
        if ($CurrentFilter.find((f) => f.name === event.data.type)) {
            return true;
        }
        return false;
    };

    const searchEvents = (event) => {
        if ($SearchQuery === '') {
            return true;
        }
        // return true if it is in any element of event.data
        if (
            Object.values(event.data).some((value) => {
                if (typeof value === 'string') {
                    return value.toLowerCase().includes($SearchQuery.toLowerCase());
                }
                return false;
            })
        ) {
            return true;
        }
        return false;
    };

    $: filteredEvents = events;

    $: futureEvents = filteredEvents
        .filter((event) => {
            const today = new Date();
            return event.data.start > today;
        })
        .sort((a, b) => {
            if (a.data.start < b.data.start) {
                return -1;
            }
            return 1;
        });

    $: pastEvents = filteredEvents
        .filter((event) => {
            const today = new Date();
            return event.data.end < today;
        })
        .sort((a, b) => {
            if (a.data.end < b.data.end) {
                return 1;
            }
            return -1;
        });

    $: if (currentEvents.length > 0) {
        EventIsOngoing.set(true);
    } else {
        EventIsOngoing.set(false);
    }

    const event_type_classes = {
        bytesize: 'success',
        hackathon: 'primary',
        talk: 'info',
        training: 'warning',
    };
    const event_type_icons = {
        bytesize: 'fa-solid fa-apple-core',
        hackathon: 'fa-solid fa-laptop-code',
        talk: 'fa-solid fa-presentation',
        training: 'fa-solid fa-chalkboard-teacher',
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

    onMount(() => {
        CurrentFilter.set(currentFilters);
        CurrentFilter.subscribe(() => {
            filteredEvents = events.filter(filterByType).filter(searchEvents);
        });

        SearchQuery.subscribe(() => {
            filteredEvents = events.filter(filterByType).filter(searchEvents);
        });
    });
</script>

<div>
    <FilterBar filter={event_types} displayStyle={[]} sortBy={[]} />
    <div class="events">
        {#if currentEvents.length > 0}
            <div class="mb-3 col-12">
                <h2><i class="fa-duotone fa-calendar-exclamation me-3" />Currently ongoing</h2>
                {#each currentEvents as event (event.id)}
                    <EventCard
                        frontmatter={event.data}
                        slug={event.slug}
                        type={event.data.type}
                        time_category="current"
                    />
                {/each}
            </div>
        {/if}
        <div>
            <h1>General events</h1>
            <div class="d-flex flex-column">
                <div class="mb-3">
                    <h2><i class="fa-duotone fa-calendar-day me-3" />Upcoming events</h2>
                    {#if futureEvents && futureEvents.length > 0}
                        {#each futureEvents as event (event.id)}
                            <EventCard
                                frontmatter={event.data}
                                slug={event.slug}
                                type={event.data.type}
                                time_category="future"
                            />
                        {/each}
                    {:else if !$SearchQuery && $CurrentFilter.length !== 0}
                        <p>Nothing in the calendar at the moment.</p>
                    {/if}
                </div>
                <div class="mb-3">
                    <h2><i class="fa-duotone fa-calendar-check me-3" />Past events</h2>
                    {#each pastEvents as event, idx (event.id)}
                        {#if hasYearChanged(pastEvents, idx)}
                            <h3 id={'year-' + event.data.start.getFullYear()}>{event.data.start.getFullYear()}</h3>
                        {/if}
                        <EventCard
                            frontmatter={event.data}
                            slug={event.slug}
                            type={event.data.type}
                            time_category="past"
                        />
                    {/each}
                </div>
            </div>
        </div>
    </div>
</div>
