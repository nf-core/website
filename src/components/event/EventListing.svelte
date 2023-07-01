<script lang="ts">
    import FilterBar from '@components/FilterBar.svelte';
    import EventCard from '@components/event/EventCard.svelte';
    import { CurrentFilter, SearchQuery, EventIsOngoing } from '@components/store';

    export let events = [];

    let filteredEvents = events;
    const filterByType = (event) => {
        if ($CurrentFilter.includes(event.data.type)) {
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

    CurrentFilter.subscribe(() => {
        filteredEvents = events.filter(filterByType).filter(searchEvents);
    });

    SearchQuery.subscribe(() => {
        filteredEvents = events.filter(filterByType).filter(searchEvents);
    });

    $: futureEvents = filteredEvents.filter((event) => {
        const today = new Date();
        return event.data.start > today;
    });

    $: pastEvents = filteredEvents
        .filter((event) => {
            const today = new Date();
            return event.data.end < today;
        })
        .reverse();
    let currentEvents: CollectionEntry<'events'>[] = [];
    $: currentEvents = filteredEvents.filter((event) => {
        const today = new Date();
        return event.data.start < today && event.data.end > today;
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
        poster: 'fa-regular fa-image',
        talk: 'fa-solid fa-presentation',
        training: 'fa-solid fa-chalkboard-teacher',
    };

    const event_types = [...new Set(events.map((event) => event.data.type))]
        .map((type) => {
            return {
                name: type,
                class: event_type_classes[type],
                icon: event_type_icons[type],
            };
        })
        .sort((a, b) => {
            if (a.name < b.name) {
                return -1;
            }
            return 1;
        });
</script>

<div>
    <FilterBar filter={event_types} displayStyle={[]} sortBy={[]} />

    {#if currentEvents.length > 0}
        <div class="mb-3 col-12 col-md-6">
            <h2><i class="fa-duotone fa-calendar-exclamation me-3" />Currently ongoing</h2>
            {#each currentEvents as event (event.id)}
                <EventCard
                    frontmatter={event.data}
                    slug={event.slug}
                    type={event.data.type}
                    time_category="current"
                    {event_type_classes}
                />
            {/each}
        </div>
    {/if}
    <div class="row">
        <div class="col-12 col-md-6">
            <h1>General events</h1>
            <div class="d-flex flex-column">
                <div class="mb-3">
                    <h2><i class="fa-duotone fa-calendar-day me-3" />Upcoming events</h2>
                    {#if futureEvents.length > 0}
                        {#each futureEvents.filter((e) => e.data.type !== 'bytesize') as event (event.id)}
                            <EventCard
                                frontmatter={event.data}
                                slug={event.slug}
                                type={event.data.type}
                                time_category="future"
                                {event_type_classes}
                            />
                        {/each}
                    {:else if $SearchQuery === '' && $CurrentFilter.length !== 0}
                        <p>Nothing in the calendar at the moment.</p>
                    {/if}
                </div>
                <div class="mb-3">
                    <h2><i class="fa-duotone fa-calendar-check me-3" />Past events</h2>
                    {#each pastEvents.filter((e) => e.data.type !== 'bytesize') as event (event.id)}
                        <EventCard
                            frontmatter={event.data}
                            slug={event.slug}
                            type={event.data.type}
                            time_category="past"
                            {event_type_classes}
                        />
                    {/each}
                </div>
            </div>
        </div>
        <div class="col-12 col-md-6">
            <h1>Bytesize talks</h1>
            <div class="d-flex flex-column">
                <div class="mb-3">
                    <h2><i class="fa-duotone fa-calendar-day me-3" />Upcoming events</h2>
                    {#if futureEvents.length > 0}
                        {#each futureEvents.filter((e) => e.data.type === 'bytesize') as event (event.id)}
                            <EventCard
                                frontmatter={event.data}
                                slug={event.slug}
                                type={event.data.type}
                                time_category="future"
                                {event_type_classes}
                            />
                        {/each}
                    {:else if $SearchQuery === '' && $CurrentFilter.length !== 0}
                        <p>No bytesize talks are scheduled at the moment.</p>
                    {/if}
                </div>
                <div class="mb-3">
                    <h2><i class="fa-duotone fa-calendar-check me-3" />Past events</h2>
                    {#each pastEvents.filter((e) => e.data.type === 'bytesize') as event (event.id)}
                        <EventCard
                            frontmatter={event.data}
                            slug={event.slug}
                            type={event.data.type}
                            time_category="past"
                            {event_type_classes}
                        />
                    {/each}
                </div>
            </div>
        </div>
    </div>
</div>
