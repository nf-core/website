<script lang="ts">
    import { ICalendar, GoogleCalendar, OutlookCalendar } from 'datebook/dist/datebook.min.mjs';
    import pkg from 'file-saver';
    const { saveAs } = pkg;

    export let frontmatter = {};
    export let add_class = 'btn-outline-success';

    let event_location = '';
    if (typeof frontmatter.locationURL === 'string') {
        event_location = frontmatter.locationURL;
    } else if (frontmatter.locationURL) {
        event_location = frontmatter.locationURL.join(', ');
    }
    const calendar_event = {
        title: frontmatter.title,
        description: frontmatter.subtitle,
        start: frontmatter.start.toISOString().includes('00:00:00')
            ? new Date(frontmatter.start.toISOString().split('T')[0])
            : frontmatter.start,
        end: frontmatter.end.toISOString().includes('00:00:00')
            ? new Date(frontmatter.end.toISOString().split('T')[0])
            : frontmatter.end,
        location: event_location,
        allDay: frontmatter.start.toISOString().includes('00:00:00'),
    };

    const removeTimeFromICSDate = (ics: string) => {
        // find DTSTART or DTEND lines and drop the time completely
        return ics
            .replace(/DTSTART:(\d{8})T\d{6}Z/g, 'DTSTART;VALUE=DATE:$1')
            .replace(/DTEND:(\d{8})T\d{6}Z/g, 'DTEND;VALUE=DATE:$1');
    };

    const removeTimeFromCalendarURL = (url: string, isGoogle: boolean = false) => {
        let params = new URL(url).searchParams;
        // Add one day to the end date for all-day events
        const endDate = new Date(calendar_event.end);
        endDate.setDate(endDate.getDate() + 1);

        params.set('allday', 'true');
        if (isGoogle) {
            // Format dates as YYYYMMDD using toLocaleDateString
            const formatDate = (date: Date) => date.toLocaleDateString('en-CA').replace(/-/g, ''); // en-CA gives YYYY-MM-DD format
            params.set('dates', formatDate(calendar_event.start) + '/' + formatDate(endDate));
        } else {
            params.set('startdt', calendar_event.start.toISOString().split('T')[0]);
            params.set('enddt', endDate.toISOString().split('T')[0]);
        }

        return url.split('?')[0] + '?' + params.toString();
    };

    let googleCalendar = new GoogleCalendar(calendar_event, true).render();
    let outlookCalendar = new OutlookCalendar(calendar_event).render();
    let ical = new ICalendar(calendar_event).render();
    if (calendar_event.allDay) {
        ical = removeTimeFromICSDate(ical);
        googleCalendar = removeTimeFromCalendarURL(googleCalendar);
        outlookCalendar = removeTimeFromCalendarURL(outlookCalendar);
    }
    const blob = new Blob([ical], { type: 'text/calendar;charset=utf-8' });
    function downloadIcal() {
        saveAs(blob, frontmatter.title.replace(/[^a-z0-9]/gi, '_').toLowerCase() + '.ics');
    }
</script>

<div class="dropdown btn-group position-relative" role="group">
    <button
        type="button"
        class={'btn dropdown-toggle ' + add_class}
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
