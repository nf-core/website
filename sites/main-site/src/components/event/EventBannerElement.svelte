<script lang="ts">
    import { formatDistanceToNow } from "date-fns";
    import ExportEventButton from "@components/event/ExportEventButton.svelte";
    import VideoButton from "@components/VideoButton.svelte";
    import { onMount } from "svelte";
    import type { CollectionEntry } from "astro:content";
    interface Props {
        events?: CollectionEntry<"events">[];
        event_time_category?: string;
        event_type_classes?: {};
        event_type_icons?: {};
    }

    let {
        events = $bindable([]),
        event_time_category = "",
        event_type_classes = {},
        event_type_icons = {},
    }: Props = $props();

    let backgroundIcon = $state("");

    const event_duration = (event) => {
        event.data.eventCountDown = formatDistanceToNow(event.data.start);
        if (event.data.startDate === event.data.endDate) {
            event.data.duration =
                new Date(event.data.start).toLocaleString("en-US", {
                    year: "numeric",
                    month: "short",
                    day: "numeric",
                    hour: "numeric",
                    minute: "numeric",
                    hour12: false,
                }) +
                "-" +
                new Date(event.data.end).toLocaleString("en-US", {
                    hour: "numeric",
                    minute: "numeric",
                    hour12: false,
                });
        } else {
            event.data.duration =
                new Date(event.data.start).toLocaleString("en-US", {
                    year: "numeric",
                    month: "short",
                    day: "numeric",
                    hour: "numeric",
                    minute: "numeric",
                    hour12: false,
                }) +
                " - " +
                new Date(event.data.end).toLocaleString("en-US", {
                    year: "numeric",
                    month: "short",
                    day: "numeric",
                    hour: "numeric",
                    minute: "numeric",
                    hour12: false,
                });
        }
    };

    events
        .map((event) => {
            if (event.data.title.toLowerCase().match("bytesize")) {
                event.data.type = "bytesize";
            }
            event_duration(event);
            return event;
        })
        .sort((a, b) => {
            return new Date(a.data.start).getTime() - new Date(b.data.start).getTime();
        });
    if (event_time_category === "upcoming") {
        backgroundIcon = "fa-alarm-clock";
        events = events
            .filter((event) => {
                let time_window = 1 * 24 * 60 * 60 * 1000;
                let event_start_unix = event.data.start.getTime();
                if (event.data.announcement?.start !== undefined) {
                    event_start_unix = new Date(event.data.announcement.start).getTime();
                    time_window = 0;
                }
                const event_end_unix = event.data.end.getTime();

                // increase time window to a week for events longer than 5 hours
                if (
                    event_end_unix - event_start_unix > 5 * 60 * 60 * 1000 &&
                    event.data.announcement?.start === undefined
                ) {
                    time_window = 7 * 24 * 60 * 60 * 1000;
                }
                if (event.data.start < new Date() && new Date() < event.data.end) {
                    // this is not an upcoming, but an ongoing event
                    return false;
                }

                if (event_start_unix < new Date().getTime() + time_window && new Date().getTime() < event_end_unix) {
                    return true;
                }
            })
            .sort((a, b) => {
                return new Date(a.data.start).getTime() - new Date(b.data.start).getTime();
            });
    } else if (event_time_category === "ongoing") {
        backgroundIcon = "fa-broadcast-tower";
        events = events
            .filter((event) => {
                return event.data.start < new Date() && new Date() < event.data.end;
            })
            .sort((a, b) => {
                return new Date(b.data.start).getTime() - new Date(a.data.start).getTime();
            });
    }

    let heading_title = $derived(
        event_time_category.charAt(0).toUpperCase() +
            event_time_category.slice(1) +
            " event" +
            (events.length > 1 ? "s" : ""),
    );
    onMount(() => {
        events.map((event) => {
            event_duration(event);
            return event;
        });
    });
</script>

{#if events.length > 0}
    <div class={event_time_category + "-event event-container border-bottom border-black-subtle"}>
        <div>
            <div class="d-none d-lg-flex">
                <div class="col-lg-4 overflow-hidden ps-3 position-relative d-flex flex-column justify-content-center">
                    <h4 class="display-4 p-2 pb-0 mb-0 flex-grow-1">{heading_title}</h4>
                    <slot name="backgroundIcon">
                        <i
                            class={`fad ${backgroundIcon} homepage-header-fa-background mt-5 ms-1 ms-xl-5`}
                            aria-hidden="true"
                        ></i>
                    </slot>
                </div>
                <div class="flex-grow-1">
                    {#each events as event (event.id)}
                        <div class="w-100 row align-items-center">
                            <div class="col-8 py-lg-2 text-lg-start">
                                <h5 class="pt-2 pb-0 pb-lg-1">
                                    <a href={"events/" + event.id + "/"} class="text-success text-decoration-none"
                                        >{event.data.title}</a
                                    >
                                    <span class="ms-1 my-auto">
                                        <span class={"badge bg-" + event_type_classes[event.data.type] + " small"}
                                            ><i class={event_type_icons[event.data.type] + " me-1"} aria-hidden="true"
                                            ></i>
                                            {event.data.type}</span
                                        >
                                    </span>
                                </h5>
                                <a href={"events/" + event.id + "/"} class="text-body text-decoration-none lead"
                                    >{@html event.data.subtitle}</a
                                >
                                {#if event.data.duration}
                                    <p class="mb-1">
                                        <a
                                            href={"events/" + event.id + "/"}
                                            class="text-secondary-emphasis text-decoration-none"
                                            >{event.data.duration}</a
                                        >
                                    </p>
                                {/if}
                            </div>

                            <div class="col-4 py-lg-2 text-start d-flex flex-column align-items-start">
                                {#if event_time_category === "upcoming"}
                                    <div class="text-nowrap ps-1">
                                        <h5>Event starts in</h5>
                                        <span class="display-6">
                                            {@html event.data.eventCountDown}
                                        </span>
                                    </div>
                                    <div class="btn-group my-2" role="group" aria-label="Event details">
                                        <a
                                            href={"events/" + event.id + "/"}
                                            class="btn btn-outline-success text-nowrap"
                                        >
                                            Event Details
                                        </a>
                                        <ExportEventButton frontmatter={event.data} />
                                    </div>
                                {/if}
                                {#if event_time_category === "ongoing"}
                                    <div class="">
                                        <div class="btn-group" role="group" aria-label="Event details">
                                            <a
                                                href={"events/" + event.id + "/"}
                                                class="btn btn-outline-success text-nowrap">Event Details</a
                                            >
                                            {#if Array.isArray(event.data?.locations) && event.data.locations.length > 0}
                                                <VideoButton urls={event.data.locations} />
                                            {/if}
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
                {#each events as event (event.id)}
                    <div class="text-center">
                        <h4 class="pt-2 pb-0">
                            <a href={"events/" + event.id + "/"} class="text-success text-decoration-none"
                                >{event.data.title}</a
                            >
                        </h4>
                        <p class="d-sm-none mb-1">
                            <a href={"events/" + event.id + "/"} class="text-body text-decoration-none"
                                >{@html event.data.subtitle}</a
                            ><span class={"badge bg-" + event_type_classes[event.data.type] + " small ms-3"}
                                ><i class={event_type_icons[event.data.type] + " me-1"} aria-hidden="true"></i>
                                {event.data.type}</span
                            >
                        </p>
                        <div class="small mb-1 mx-3 d-flex flex-column">
                            <a
                                href={"events/" + event.id + "/"}
                                class="text-secondary-emphasis text-decoration-none mb-2">{event.data.duration}</a
                            >
                            <div class="btn-group text-nowrap" role="group" aria-label="Event details">
                                <a href={"events/" + event.id + "/"} class="btn btn-outline-success"> Event Details </a>
                                {#if event_time_category === "upcoming"}
                                    <ExportEventButton frontmatter={event.data} />
                                {/if}
                                {#if event_time_category === "ongoing" && Array.isArray(event.data?.locations) && event.data.locations.length > 0}
                                    <VideoButton urls={event.data.locations} />
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
    hr:last-child {
        display: none;
    }
</style>
