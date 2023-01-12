<script>
    import ExportEventButton from './ExportEventButton.svelte';
    import VideoButton from './VideoButton.svelte';
    export let events = [];
    export let event_type_classes = {};
    export let event_type_icons = {};

    const upcomingEvents = events.filter((event) => {
        let time_window = 24 * 60 * 60 * 1000;
        const event_start_unix = event.data.start.getTime();
        const event_end_unix = event.data.end.getTime();

        // increase time window to a week for events longer than 5 hours
        if (event_end_unix - event_start_unix > 5 * 60 * 60 * 1000) {
            time_window = 7 * 24 * 60 * 60 * 1000;
        }
        if (event.data.start < new Date() && new Date() < event.data.end) {
            return false;
        }

        // if (event_start_unix < new Date().getTime() + time_window && new Date().getTime() < event_end_unix) {
        // TODO: uncoment above line and remove below line, which is only for testing purposes
        if (new Date().getTime() < event_end_unix) {
            return true;
        }
    });

    const heading_title = upcomingEvents.length > 1 ? 'Upcoming events' : 'Upcoming event';

    // countdown function to event start which updates every second, by making `now` and the function reactive
    let now = '';
    setInterval(() => (now = new Date().getTime()), 1000);
    $: countdown = (event_start, short = false) => {
        let timeLeftString = '';
        const distance = event_start.getTime() - now;
        // Time calculations for days, hours, minutes and seconds
        let days = Math.floor(distance / (1000 * 60 * 60 * 24));
        const hours = Math.floor((distance % (1000 * 60 * 60 * 24)) / (1000 * 60 * 60));
        const minutes = Math.floor((distance % (1000 * 60 * 60)) / (1000 * 60));
        const seconds = Math.floor((distance % (1000 * 60)) / 1000);
        // Display the result in the element with id="demo"
        timeLeftString = days + ' days,\n ' + hours + 'h ' + minutes + 'm';
        if (short) {
            timeLeftString = days + 'd ' + hours + 'h ' + minutes + 'm';
        }
        // If the count down is finished, write some text
        if (days < 0) {
            timeLeftString = hours + 'h ' + minutes + 'm ' + seconds + 's';
        }
        if (distance < 0) {
            timeLeftString = 'Event started';
        }
        return timeLeftString;
    };
</script>

{#if upcomingEvents.length > 0}
    <div class="upcoming-event">
        <div>
            <div class="d-none d-lg-flex">
                <div class="col-lg-4 overflow-hidden">
                    <h4 class="display-4 pt-2">{heading_title}</h4>
                    <i class="fad fa-alarm-clock homepage-header-fa-background ms-1 ms-xl-5" aria-hidden="true" />
                </div>
                <div>
                    {#each upcomingEvents as event}
                        <div class="row d-flex align-items-center">
                            <div class="col pt-lg-3 pb-lg-3 text-lg-start border-bottom border-black-subtle">
                                <h5 class="pt-2 pb-0 pb-lg-1">
                                    <a href={'events/' + event.slug} class="text-success text-decoration-none"
                                        >{event.data.title}</a
                                    >
                                </h5>
                                <p class="lead">
                                    <a href={'events/' + event.slug} class="text-body text-decoration-none"
                                        >{event.data.subtitle}</a
                                    >
                                </p>

                                <p class="">
                                    <a href={'events/' + event.slug} class="text-secondary text-decoration-none"
                                        >{event.data.duration}</a
                                    ><span class="ms-1">
                                        <span class={'badge bg-' + event_type_classes[event.data.type] + ' small'}
                                            ><i
                                                class={event_type_icons[event.data.type] + ' me-1'}
                                                aria-hidden="true"
                                            />
                                            {event.data.type}</span
                                        >
                                    </span>
                                </p>
                                <div class="btn-group" role="group" aria-label="Event details">
                                    <a href={'events/' + event.slug} class="btn btn-outline-success text-nowrap">
                                        Event Details
                                    </a>
                                    <ExportEventButton frontmatter={event.data} />
                                </div>
                            </div>

                            <div class="col-lg-4 col-xl-3">
                                <div class="text-center">
                                    <h5>Event starts in:</h5>
                                    <span class="display-5">
                                        {countdown(event.data.start)}
                                    </span>
                                </div>
                            </div>
                        </div>
                    {/each}
                </div>
            </div>
            <div class="d-lg-none">
                <div class="pt-2 pb-1 mb-2 overflow-hidden mainpage-subheader-heading-header bg-body-tertiary">
                    <h5 class="pt-2 font-weight-light text-center text-sucess">{heading_title}</h5>
                </div>
                {#each upcomingEvents as event}
                    <div class="text-center border-bottom border-black-subtle">
                        <h4 class="pt-2 pb-0">
                            <a href={'events/' + event.slug} class="text-success text-decoration-none"
                                >{event.data.title}</a
                            >
                        </h4>
                        <p class="d-sm-none mb-2">
                            <a href={'events/' + event.slug} class="text-body text-decoration-none"
                                >{event.data.subtitle}</a
                            ><span class={'badge bg-' + event_type_classes[event.data.type] + ' small ms-3'}
                                            ><i
                                                class={event_type_icons[event.data.type] + ' me-1'}
                                                aria-hidden="true"
                                            />
                                            {event.data.type}</span
                                        >
                        </p>
                        <div class="small mb-2">
                            <a href={'events/' + event.slug} class="text-secondary text-decoration-none"
                                >{event.data.duration}</a
                            >
                            <div class="btn-group" role="group" aria-label="Event details">
                                <a href={'events/' + event.slug} class="btn btn-outline-success text-nowrap">
                                    Event Details
                                </a>
                                <ExportEventButton frontmatter={event.data} />
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
        margin-top: -2.5rem;
        font-size: 16em;
        opacity: 0.2;
    }
    .row:last-child .col {
        border-bottom: 0 !important;
        padding-bottom: 1pt !important;
    }
</style>
