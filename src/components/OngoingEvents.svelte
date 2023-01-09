<script>
    import { EventIsOngoing } from './store.js';
    import VideoButton from './VideoButton.svelte';

    export let events = [];
    export let event_type_classes = {};
    export let event_type_icons = {};
    const ongoingEvents = events
        .filter((event) => {
            if (event.data.start < new Date() && new Date() < event.data.end) {
                EventIsOngoing.set(true);
                return true;
            }
        })
        .sort((a, b) => {
            return a.data.start - b.data.start;
        });
    console.log('ongoing', ongoingEvents[0].data);
</script>

{#if ongoingEvents.length > 0}
    <div class="ongoing-event">
        <div>
            <div class="d-none d-lg-flex">
                <div class="col-lg-3 overflow-hidden">
                    <h4 class="display-4 pt-2">Ongoing event</h4>
                    <i class="fad fa-broadcast-tower homepage-header-fa-background ms-1 ms-xl-5" aria-hidden="true" />
                </div>
                <div>
                    {#each ongoingEvents as event}
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
                            </div>

                            <div class="col-lg-4 col-xl-3">
                                <div class="pt-lg-5">
                                    <VideoButton urls={event.data.location_url} />
                                    <a href={'events/' + event.data.slug} class="btn btn-outline-success mb-2"
                                        >Event Details</a
                                    >
                                </div>
                            </div>
                        </div>
                    {/each}
                </div>
            </div>
            <div class="d-lg-none">
                <div class="pt-2 pb-1 mb-2 overflow-hidden mainpage-subheader-heading-header">
                    <h5 class="pt-2 font-weight-light text-center">Ongoing event</h5>
                </div>
                {#each ongoingEvents as event}
                    <div class="row">
                        <div class="col text-center border-bottom border-black-subtle">
                            <h4 class="pt-2 pb-0">
                                <a href={'events/' + event.data.slug} class="text-success text-decoration-none"
                                    >{event.data.title}</a
                                >
                            </h4>
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
    @import '../styles/_variables.scss';
    .mainpage-subheader-heading {
        background-color: $gray-200;
        margin-top: -3.5rem;
        padding: 4rem 0 1rem;
    }
    @include media-breakpoint-down(md) {
        .mainpage-subheader-heading {
            padding-top: 2rem;
        }
    }

    .mainpage-subheader-heading-header {
        background-color: $gray-300;
        box-shadow: 0 0 10px rgba($black, 0.2);
    }
    .mainpage-subheader-heading a.subheader-link {
        text-decoration: none;
        color: $gray-800;
    }
    .homepage-header-fa-background {
        position: absolute;
        margin-top: -2.5rem;
        font-size: 16em;
        opacity: 0.2;
    }
    .row:last-child .col {
        border-bottom: 0;
    }
</style>
