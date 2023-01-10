<script>
    import ExportEventButton from "./ExportEventButton.svelte";
import VideoButton from "./VideoButton.svelte";
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
        if (new Date().getTime() < event_end_unix ) {
            return true;
        }
    });

    const countdown = (event_start) => {
        let timeLeftString =''
        setInterval((timeLeftString) => {
        const now = new Date();
        const timeLeft = event_start - now;
        const days = Math.floor(timeLeft / (1000 * 60 * 60 * 24));
        const hours = Math.floor((timeLeft % (1000 * 60 * 60 * 24)) / (1000 * 60 * 60));
        const minutes = Math.floor((timeLeft % (1000 * 60 * 60)) / (1000 * 60));
        const seconds = Math.floor((timeLeft % (1000 * 60)) / 1000);
        timeLeftString = `${days} days, ${hours}h:${minutes}m`;
        if (days <= 0){
            timeLeftString = `${hours}h:${minutes}m:${seconds}s`;
        }
        return timeLeftString;
    }, 1000);
    console.log('timeLeftString',timeLeftString);
}

console.log('countdown',countdown(upcomingEvents[0].data.start));
</script>

{#if upcomingEvents.length > 0}
    <div class="upcoming-event">
        <div>
            <div class="d-none d-lg-flex">
                <div class="col-lg-3 overflow-hidden">
                    <h4 class="display-4 pt-2">Upcoming event</h4>
                    <i class="fad fa-alarm-clock homepage-header-fa-background ms-1 ms-xl-5" aria-hidden="true" />
                </div>
                <div>
                    {#each upcomingEvents as event}
                        <div class="row">
                            <div class="col pt-lg-3 pb-lg-3 text-lg-start border-bottom border-black-subtle">
                                <h5 class="pt-2 pb-0 pb-lg-1">
                                    <a href={'events/' + event.data.slug} class="text-success text-decoration-none"
                                        >{event.data.title}</a
                                    >
                                </h5>
                                <p class="lead">
                                    <a href={'events/' + event.data.slug} class="text-body text-decoration-none"
                                        >{event.data.subtitle}</a
                                    >
                                </p>

                                <p class="">
                                    <a href={'events/' + event.data.slug} class="text-secondary text-decoration-none"
                                        >{event.data.duration}</a
                                    ><span class="ms-1">
                                        <span class={'badge bg-' + event_type_classes[event.data.type] + ' small'}
                                            ><i
                                                class={event_type_icons[event.data.type] + ' me-1'}
                                                aria-hidden="true"
                                            />
                                            Talk</span
                                        >
                                    </span>
                                </p>
                                    <div class="btn-group" role="group" aria-label="Event details">
                                        <button  type="button" href={'events/' + event.data.slug} class="btn btn-outline-success mb-2"
                                        >
                                        Event Details
                                    </button>
                                        <ExportEventButton frontmatter={event.data} />
                                    </div>
                            </div>

                            <div class="col-lg-4 col-xl-3">
                                <div class="pt-lg-5">
                                    <h5>
                                        Event countdown:
                                    </h5>
                                    <p class="lead">
                                        {countdown(event.data.start)}
                                </div>
                            </div>
                        </div>
                    {/each}
                </div>
            </div>
            <div class="d-lg-none">
                <div class="pt-2 pb-1 mb-2 overflow-hidden mainpage-subheader-heading-header">
                    <h5 class="pt-2 font-weight-light text-center">Upcoming event</h5>
                </div>
                {#each upcomingEvents as event}
                    <div class="row">
                        <div class="col text-center border-bottom border-black-subtle">
                            <h6 class="pt-2 pb-0">
                                <a href={'events/' + event.data.slug} class="text-success text-decoration-none"
                                    >{event.data.title}</a
                                >
                            </h6>
                            <p class="d-sm-none mb-2">
                                <a href={'events/' + event.data.slug} class="text-body text-decoration-none"
                                    >{event.data.subtitle}</a
                                >
                            </p>
                            <p class="small mb-2">
                                <a href={'events/' + event.data.slug} class="text-secondary text-decoration-none"
                                    >{event.data.duration}</a
                                ><VideoButton urls={event.data.location_url} />
                                <a href={'events/' + event.data.slug} class="btn btn-outline-success mb-2"
                                    >Event Details</a
                                >
                            </p>
                        </div>
                    </div>
                {/each}
            </div>
        </div>
    </div>
{/if}

<style lang="scss">
</style>
