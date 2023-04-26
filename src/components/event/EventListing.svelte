<script>
    export let events = [];

    import EventCard from '@components/event/EventCard.svelte';
    import FilterBar from '@components/FilterBar.svelte';
    import { CurrentFilter, SearchQuery, EventIsOngoing } from '@components/store.js';

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

    $: filteredEvents = events.filter(filterByType).filter(searchEvents);

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
    let currentEvents = [];
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
        poster: 'danger',
        talk: 'success',
        tutorial: 'info',
        training: 'warning',
    };
    const event_type_icons = {
        bytesize: 'fa-solid fa-apple-core',
        hackathon: 'fa-solid fa-laptop-code',
        poster: 'fa-regular fa-image',
        talk: 'fa-solid fa-presentation',
        tutorial: 'fa-solid fa-graduation-cap',
        training: 'fa-solid fa-chalkboard-teacher',
    };

    const event_types = [...new Set(events.map((event) => event.data.type))].map((type) => {
        return {
            name: type,
            class: event_type_classes[type],
            icon: event_type_icons[type],
        };
    });
</script>

<div>
    <FilterBar filter={event_types} displayStyle={[]} sortBy={[]} />
    <div>
        {#if currentEvents.length > 0}
            <div class="mb-3 mb-md-5 col">
                <h2><i class="fa-duotone fa-calendar-exclamation me-3" />Currently ongoing</h2>
                {#each currentEvents as event (event.id)}
                    <EventCard
                        frontmatter={event.data}
                        slug={event.slug.startsWith('/events/') ? event.slug.replace('/events', '') : event.slug}
                        type_class={event_type_classes[event.data.type]}
                        time_category="current"
                    />
                {/each}
            </div>
        {/if}
        <div class="row">
            <div class="mb-3 mb-md-5 col-12 col-md-6">
                <h2><i class="fa-duotone fa-calendar-day me-3" />Upcoming events</h2>
                {#if futureEvents.length > 0}
                    {#each futureEvents as event (event.id)}
                        <EventCard
                            frontmatter={event.data}
                            slug={event.slug}
                            type_class={event_type_classes[event.data.type]}
                            time_category="future"
                        />
                    {/each}
                {:else if $SearchQuery === '' && $CurrentFilter.length !== 0}
                    <p>No upcoming events at the moment</p>
                {/if}
            </div>
            <div class="mb-3 mb-md-5 col-12 col-md-6">
                <h2><i class="fa-duotone fa-calendar-check me-3" />Past events</h2>
                {#each pastEvents as event (event.id)}
                    <EventCard
                        frontmatter={event.data}
                        slug={event.slug}
                        type_class={event_type_classes[event.data.type]}
                        time_category="past"
                    />
                {/each}
            </div>
        </div>
    </div>
</div>
