<script>
    import { ICalendar, GoogleCalendar, OutlookCalendar } from 'datebook/dist/datebook.min.mjs';
    import { saveAs } from 'file-saver';

    export let frontmatter = {};
    export let add_class="btn-outline-success";

    let event_location = '';
    if (typeof frontmatter.location_url === 'string') {
        event_location = frontmatter.location_url;
    } else {
        event_location = frontmatter.location_url.join(', ');
    }
    const calendar_event = {
        title: frontmatter.title,
        description: frontmatter.subtitle,
        start: frontmatter.start,
        end: frontmatter.end,
        location: event_location,
    };

    const googleCalendar = new GoogleCalendar(calendar_event).render();
    const outlookCalendar = new OutlookCalendar(calendar_event).render();
    let ical = new ICalendar(calendar_event).render();
    const blob = new Blob([ical], { type: 'text/calendar;charset=utf-8' });
    function downloadIcal() {
        saveAs(blob, frontmatter.title.replace(/[^a-z0-9]/gi, '_').toLowerCase() + '.ics');
    }
</script>

<div class="dropwdown btn-group" role="group">
    <button
        type="button"
        class={"btn dropdown-toggle "+add_class}
        href="#"
        data-bs-toggle="dropdown"
        aria-haspopup="true"
        aria-expanded="false"
    >
        <i class="far fa-calendar-plus me-1" aria-hidden="true" /> Export event
    </button>
    <div class="dropdown-menu">
        <button class="dropdown-item" on:click={() => downloadIcal()}> Download iCal Event</button>
        <a class="dropdown-item" href={googleCalendar} target="_blank" rel="noreferrer"> Add to Google Calendar</a>
        <a class="dropdown-item" href={outlookCalendar} target="_blank" rel="noreferrer"> Add to Microsoft Outlook</a>
    </div>
</div>
