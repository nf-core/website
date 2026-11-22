import type { CollectionEntry } from "astro:content";

type EventData = CollectionEntry<"events">["data"];
export type EventType = EventData["type"];

export const eventTypeClasses: Record<EventType, string> = {
    bytesize: "success",
    hackathon: "primary",
    talk: "info",
    training: "warning",
};

export const eventTypeIcons: Record<EventType, string> = {
    bytesize: "fa-solid fa-apple-core",
    hackathon: "fa-solid fa-laptop-code",
    talk: "fa-solid fa-presentation",
    training: "fa-solid fa-chalkboard-teacher",
};

const HOUR = 60 * 60 * 1000;
const DAY = 24 * HOUR;

/**
 * Whether an event should be shown in the homepage event banner: it hasn't ended yet and it starts within its
 * announcement window. By default the window is 1 day before the start, or 7 days for events longer than 5 hours.
 * An `announcement.start` without `announcement.text` replaces that window (announcements with text go to the
 * navbar AnnouncementBanner instead).
 * `offSet` widens the window, e.g. for the build-time prefilter so events don't miss the daily static build.
 */
export function isInBannerWindow(data: EventData, now = Date.now(), offSet = 0): boolean {
    if (!data.start || !data.end) return false;
    const end = data.end.getTime();

    let start: number;
    let window: number;
    if (data.announcement?.start && !data.announcement.text) {
        start = new Date(data.announcement.start).getTime();
        window = 0;
    } else {
        start = data.start.getTime();
        window = end - start > 5 * HOUR ? 7 * DAY : DAY;
    }

    return start < now + window + offSet && now < end;
}
