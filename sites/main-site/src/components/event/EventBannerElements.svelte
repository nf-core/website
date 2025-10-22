<script lang="ts">
    import { formatDistanceToNow } from "date-fns";
    import { onMount } from "svelte";

    type EventState = { element: HTMLElement; visible: boolean; countdown?: string };

    let ongoingContainer: HTMLElement | undefined = $state();
    let upcomingContainer: HTMLElement | undefined = $state();

    let ongoingEvents: EventState[] = $state([]);
    let upcomingEvents: EventState[] = $state([]);
    let showBanner = $state(false);

    let ongoingCount = $derived(ongoingEvents.filter((e) => e.visible).length);
    let upcomingCount = $derived(upcomingEvents.filter((e) => e.visible).length);
    let totalCount = $derived(ongoingCount + upcomingCount);

    let ongoingHeading = $derived("Ongoing event" + (ongoingCount !== 1 ? "s" : ""));
    let upcomingHeading = $derived("Upcoming event" + (upcomingCount !== 1 ? "s" : ""));

    function isEventOngoing(start: Date, end: Date, now: Date): boolean {
        return start < now && now < end;
    }

    function isEventUpcoming(
        start: Date,
        end: Date,
        now: Date,
        announcementStart?: string,
    ): { visible: boolean; countdown?: string } {
        let event_start_unix = start.getTime();
        let time_window = 1 * 24 * 60 * 60 * 1000;

        if (announcementStart) {
            event_start_unix = new Date(announcementStart).getTime();
            time_window = 0;
        }

        const event_end_unix = end.getTime();

        if (event_end_unix - event_start_unix > 5 * 60 * 60 * 1000 && !announcementStart) {
            time_window = 7 * 24 * 60 * 60 * 1000;
        }

        if (isEventOngoing(start, end, now)) {
            return { visible: false };
        }

        const visible = event_start_unix < now.getTime() + time_window && now.getTime() < event_end_unix;

        return {
            visible,
            countdown: visible ? formatDistanceToNow(start) : undefined,
        };
    }

    function processEvents() {
        const now = new Date();

        if (ongoingContainer) {
            const elements = Array.from(ongoingContainer.querySelectorAll<HTMLElement>("[data-event-item]"));
            ongoingEvents = elements.map((element) => {
                const start = new Date(element.dataset.eventStart!);
                const end = new Date(element.dataset.eventEnd!);
                return {
                    element,
                    visible: isEventOngoing(start, end, now),
                };
            });
        }

        if (upcomingContainer) {
            const elements = Array.from(upcomingContainer.querySelectorAll<HTMLElement>("[data-event-item]"));
            upcomingEvents = elements.map((element) => {
                const start = new Date(element.dataset.eventStart!);
                const end = new Date(element.dataset.eventEnd!);
                const announcementStart = element.dataset.eventAnnouncementStart;

                const { visible, countdown } = isEventUpcoming(start, end, now, announcementStart);

                return { element, visible, countdown };
            });
        }
    }

    function updateEventVisibility(events: EventState[]) {
        events.forEach(({ element, visible, countdown }) => {
            element.style.display = visible ? "" : "none";
            if (countdown) {
                const countdownEl = element.querySelector("[data-event-countdown]");
                if (countdownEl) {
                    countdownEl.textContent = countdown;
                }
            }
        });
    }

    function updateHeading(container: HTMLElement | undefined, heading: string) {
        if (!container) return;

        const headingEl = container.querySelector("[data-event-heading]");
        if (headingEl) {
            headingEl.textContent = heading;
        }
    }

    onMount(() => {
        processEvents();
    });

    $effect(() => {
        updateEventVisibility(ongoingEvents);
        updateHeading(ongoingContainer, ongoingHeading);
    });

    $effect(() => {
        updateEventVisibility(upcomingEvents);
        updateHeading(upcomingContainer, upcomingHeading);
    });

    $effect(() => {
        // Trigger animation when any events are visible
        if (totalCount > 0 && !showBanner) {
            showBanner = true;
        }
    });
</script>

<div class="banner-wrapper" class:visible={showBanner} style:--event-count={totalCount}>
    <div bind:this={ongoingContainer} data-event-category="ongoing" class:d-none={ongoingCount === 0}>
        <slot name="ongoing" />
    </div>

    <div bind:this={upcomingContainer} data-event-category="upcoming" class:d-none={upcomingCount === 0}>
        <slot name="upcoming" />
    </div>
</div>

<style>
    .banner-wrapper {
        max-height: 0;
        overflow: hidden;
        /*transition: max-height 0.3s ease-in-out;*/
    }

    .banner-wrapper.visible {
        /* Desktop: 75px per event */
        max-height: calc(var(--event-count) * 75px);
    }

    /* Mobile: 200px per event (stacked layout) */
    /* Bootstrap lg breakpoint (992px) minus 0.02px */
    @media (max-width: 991.98px) {
        .banner-wrapper.visible {
            max-height: calc(var(--event-count) * 200px);
        }
    }
</style>
