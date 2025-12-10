<script lang="ts">
    import { ICalendar, GoogleCalendar, OutlookCalendar } from "datebook/dist/datebook.min.mjs";
    import pkg from "file-saver";
    const { saveAs } = pkg;

    interface Props {
        frontmatter?: any;
        add_class?: string;
    }

    let { frontmatter = {}, add_class = "btn-outline-success" }: Props = $props();

    const event_location = $derived(
        typeof frontmatter.locationURL === "string"
            ? frontmatter.locationURL
            : frontmatter.locationURL
              ? frontmatter.locationURL.join(", ")
              : "",
    );

    const calendar_event = $derived({
        title: frontmatter.title,
        description: frontmatter.subtitle,
        start: frontmatter.start.toISOString().includes("00:00:00")
            ? new Date(frontmatter.start.toISOString().split("T")[0])
            : frontmatter.start,
        end: frontmatter.end.toISOString().includes("00:00:00")
            ? new Date(frontmatter.end.toISOString().split("T")[0])
            : frontmatter.end,
        location: event_location,
        allDay: frontmatter.start.toISOString().includes("00:00:00"),
    });

    const removeTimeFromICSDate = (ics: string) => {
        // find DTSTART or DTEND lines and drop the time completely
        return ics
            .replace(/DTSTART:(\d{8})T\d{6}Z/g, "DTSTART;VALUE=DATE:$1")
            .replace(/DTEND:(\d{8})T\d{6}Z/g, "DTEND;VALUE=DATE:$1");
    };

    const removeTimeFromCalendarURL = (url: string, isGoogle: boolean = false) => {
        let params = new URL(url).searchParams;
        // Add one day to the end date for all-day events
        const endDate = new Date(calendar_event.end);
        endDate.setDate(endDate.getDate() + 1);

        params.set("allday", "true");
        if (isGoogle) {
            // Format dates as YYYYMMDD using toLocaleDateString
            const formatDate = (date: Date) => date.toLocaleDateString("en-CA").replace(/-/g, ""); // en-CA gives YYYY-MM-DD format
            params.set("dates", formatDate(calendar_event.start) + "/" + formatDate(endDate));
        } else {
            params.set("startdt", calendar_event.start.toISOString().split("T")[0]);
            params.set("enddt", endDate.toISOString().split("T")[0]);
        }

        return url.split("?")[0] + "?" + params.toString();
    };

    // Initialize the calendar objects with derived values
    const calendarData = $derived.by(() => {
        const google = new GoogleCalendar(calendar_event, true).render();
        const outlook = new OutlookCalendar(calendar_event).render();
        let ical = new ICalendar(calendar_event).render();

        // Apply modifications if needed
        if (calendar_event.allDay) {
            ical = removeTimeFromICSDate(ical);
            return {
                google: removeTimeFromCalendarURL(google, true),
                outlook: removeTimeFromCalendarURL(outlook),
                ical: ical,
            };
        }

        return { google, outlook, ical };
    });

    // Create the blob from the ical data
    const icalBlob = $derived(new Blob([calendarData.ical], { type: "text/calendar;charset=utf-8" }));

    function downloadIcal() {
        saveAs(icalBlob, frontmatter.title.replace(/[^a-z0-9]/gi, "_").toLowerCase() + ".ics");
    }
</script>

<div class="dropdown btn-group position-relative" role="group">
    <button
        type="button"
        class={"btn dropdown-toggle " + add_class}
        data-bs-toggle="dropdown"
        aria-haspopup="true"
        aria-expanded="false"
    >
        <i class="far fa-calendar-plus me-1" aria-hidden="true"></i> Export event
    </button>
    <div class="dropdown-menu">
        <button class="dropdown-item" onclick={() => downloadIcal()}> Download iCal Event</button>
        <a class="dropdown-item" href={calendarData.google} target="_blank" rel="noreferrer"> Add to Google Calendar</a>
        <a class="dropdown-item" href={calendarData.outlook} target="_blank" rel="noreferrer">
            Add to Microsoft Outlook</a
        >
    </div>
</div>
