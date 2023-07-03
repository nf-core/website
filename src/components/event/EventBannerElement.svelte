<script lang="ts">
    import { EventIsOngoing } from '@components/store';
    import { formatDistanceToNow } from 'date-fns';
    import ExportEventButton from '@components/event/ExportEventButton.svelte';
    import VideoButton from '@components/VideoButton.svelte';
    export let events = [];
    export let event_time_category: string = '';

    export let event_type_classes: {}[] = [{}];
    export let event_type_icons: {}[] = [{}];

    let now = new Date().getTime();

    let backgroundIcon = '';
    console.log('events before processing', events);
    events
        .map((event) => {
            if (event.data.title.toLowerCase().match('bytesize')) {
                event.data.type = 'bytesize';
            }
            if (event.data.start_date === event.data.end_date) {
                event.data.duration =
                    event.data.start.toLocaleString('en-US', {
                        year: 'numeric',
                        month: 'short',
                        day: 'numeric',
                        hour: 'numeric',
                        minute: 'numeric',
                        hour12: false,
                    }) +
                    '-' +
                    event.data.end.toLocaleString('en-US', {
                        hour: 'numeric',
                        minute: 'numeric',
                        hour12: false,
                    });
            } else {
                event.data.duration =
                    event.data.start.toLocaleString('en-US', {
                        year: 'numeric',
                        month: 'short',
                        day: 'numeric',
                        hour: 'numeric',
                        minute: 'numeric',
                        hour12: false,
                    }) +
                    ' - ' +
                    event.data.end.toLocaleString('en-US', {
                        year: 'numeric',
                        month: 'short',
                        day: 'numeric',
                        hour: 'numeric',
                        minute: 'numeric',
                        hour12: false,
                    });
            }
        })
        .sort((a, b) => {
            return new Date(a.data.start) - new Date(b.data.start);
        });
    console.log('events after processing', events);
    if (event_time_category === 'upcoming') {
        backgroundIcon = 'fa-alarm-clock';
        events = events.filter((event) => {
            let time_window = 1 * 24 * 60 * 60 * 1000; //TODO: change back to 1 day
            const event_start_unix = event.data.start.getTime();
            const event_end_unix = event.data.end.getTime();

            // increase time window to a week for events longer than 5 hours
            if (event_end_unix - event_start_unix > 5 * 60 * 60 * 1000) {
                time_window = 7 * 24 * 60 * 60 * 1000; //TODO: change back to 7 days
            }
            if (event.data.start < new Date() && new Date() < event.data.end) {
                return false;
            }

            if (event_start_unix < now + time_window && now < event_end_unix) {
                return true;
            }
        });
    } else if (event_time_category === 'ongoing') {
        backgroundIcon = 'fa-broadcast-tower';
        events = events
            .filter((event) => {
                return event.data.start < new Date() && new Date() < event.data.end;
            })
            .sort((a, b) => {
                return a.data.start - b.data.start;
            });

        if (events.length > 0) {
            EventIsOngoing.set(true);
        } else {
            EventIsOngoing.set(false);
        }
    }

    let heading_title = event_time_category.charAt(0).toUpperCase() + event_time_category.slice(1) + ' event';
    heading_title = events.length > 1 ? heading_title + 's' : heading_title;

    // countdown function to event start which updates every second, by making `now` and the function reactive

    setInterval(() => (now = new Date().getTime()), 1000);
    $: countdown = (event_start) => {
        const timeLeftString = formatDistanceToNow(event_start);
        return timeLeftString;
    };
</script>

{#if events.length > 0}
    <div class={event_time_category + '-event event-container border-bottom border-black-subtle'}>
        <div>
            <div class="d-none d-lg-flex">
                <div class="col-lg-4 overflow-hidden ps-3 position-relative d-flex flex-column justify-content-center">
                    <h4 class="display-4 p-2 flex-grow-1">{heading_title}</h4>
                    <i
                        class={`fad ${backgroundIcon} homepage-header-fa-background mt-5 ms-1 ms-xl-5`}
                        aria-hidden="true"
                    />
                </div>
                <div class="flex-grow-1">
                    {#each events as event}
                        <div class="w-100 row align-items-center">
                            <div class="col-9 pt-lg-3 pb-lg-3 text-lg-start">
                                <h5 class="pt-2 pb-0 pb-lg-1">
                                    <a href={'events/' + event.slug + '/'} class="text-success text-decoration-none"
                                        >{event.data.title}</a
                                    >
                                    <span class="ms-1 my-auto">
                                        <span class={'badge bg-' + event_type_classes[event.data.type] + ' small'}
                                            ><i
                                                class={event_type_icons[event.data.type] + ' me-1'}
                                                aria-hidden="true"
                                            />
                                            {event.data.type}</span
                                        >
                                    </span>
                                </h5>
                                <p class="lead">
                                    <a href={'events/' + event.slug + '/'} class="text-body text-decoration-none"
                                        >{event.data.subtitle}</a
                                    >
                                </p>
                                {#if event.data.duration}
                                    <p class="">
                                        <a
                                            href={'events/' + event.slug + '/'}
                                            class="text-secondary-emphasis text-decoration-none"
                                            >{event.data.duration}</a
                                        >
                                    </p>
                                {/if}
                                {#if event_time_category === 'upcoming'}
                                    <div class="btn-group" role="group" aria-label="Event details">
                                        <a
                                            href={'events/' + event.slug + '/'}
                                            class="btn btn-outline-success text-nowrap"
                                        >
                                            Event Details
                                        </a>
                                        <ExportEventButton frontmatter={event.data} />
                                    </div>
                                {/if}
                            </div>

                            <div class="col-3">
                                {#if event_time_category === 'upcoming'}
                                    <div class="text-nowrap">
                                        <h5>Event starts in</h5>
                                        <span class="display-6">
                                            {@html countdown(event.data.start)}
                                        </span>
                                    </div>
                                {/if}
                                {#if event_time_category === 'ongoing'}
                                    <div class="">
                                        <div class="btn-group" role="group" aria-label="Event details">
                                            <a
                                                href={'events/' + event.slug + '/'}
                                                class="btn btn-outline-success text-nowrap">Event Details</a
                                            >
                                            <VideoButton urls={event.data.location_url} />
                                        </div>
                                    </div>
                                {/if}
                            </div>
                        </div>
                        <hr class="mx-4 my-0 py-0" />
                    {/each}
                </div>
            </div>
            <div class="d-lg-none">
                <div class="pt-2 pb-1 mb-2 overflow-hidden mainpage-subheader-heading-header bg-body-tertiary">
                    <h5 class="pt-2 font-weight-light text-center text-sucess">{heading_title}</h5>
                </div>
                {#each events as event}
                    <div class="text-center">
                        <h4 class="pt-2 pb-0">
                            <a href={'events/' + event.slug + '/'} class="text-success text-decoration-none"
                                >{event.data.title}</a
                            >
                        </h4>
                        <p class="d-sm-none mb-1">
                            <a href={'events/' + event.slug + '/'} class="text-body text-decoration-none"
                                >{event.data.subtitle}</a
                            ><span class={'badge bg-' + event_type_classes[event.data.type] + ' small ms-3'}
                                ><i class={event_type_icons[event.data.type] + ' me-1'} aria-hidden="true" />
                                {event.data.type}</span
                            >
                        </p>
                        <div class="small mb-1 mx-3 d-flex flex-column">
                            <a
                                href={'events/' + event.slug + '/'}
                                class="text-secondary-emphasis text-decoration-none mb-2">{event.data.duration}</a
                            >
                            <div class="btn-group text-nowrap" role="group" aria-label="Event details">
                                <a href={'events/' + event.slug + '/'} class="btn btn-outline-success">
                                    Event Details
                                </a>
                                {#if event_time_category === 'upcoming'}
                                    <ExportEventButton frontmatter={event.data} />
                                {/if}
                                {#if event_time_category === 'ongoing'}
                                    <VideoButton urls={event.data.location_url} />
                                {/if}
                            </div>
                        </div>
                    </div>
                {/each}
            </div>
        </div>
    </div>
{/if}

<style lang="scss">
    .homepage-header-fa-background {
        position: absolute;
        font-size: 14em;
        opacity: 0.2;
    }
    hr:last-child {
        display: none;
    }
</style>
