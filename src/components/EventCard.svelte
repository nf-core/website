<script>
    import { ICalendar, GoogleCalendar, OutlookCalendar } from 'datebook/dist/datebook.min.mjs'
    import { saveAs } from 'file-saver';
    export let frontmatter = {
        title: '',
        subtitle: '',
        start: '',
        start_date: '',
        end: '',
        end_date: '',
        type: '',
    };
    export let slug = '';
    export let type_class = '';
    export let time_category = '';

    let event_date;
    if (frontmatter.start_date === frontmatter.end_date) {
        event_date =
            frontmatter.start.toLocaleString('en-US', {
                year: 'numeric',
                month: 'short',
                day: 'numeric',
                hour: 'numeric',
                minute: 'numeric',
                hour12: false,
            }) +
            '-' +
            frontmatter.end.toLocaleString('en-US', {
                hour: 'numeric',
                minute: 'numeric',
                hour12: false,
            });
    } else {
        event_date =
            frontmatter.start.toLocaleString('en-US', {
                year: 'numeric',
                month: 'short',
                day: 'numeric',
                hour: 'numeric',
                minute: 'numeric',
                hour12: false,
            }) +
            ' - ' +
            frontmatter.end.toLocaleString('en-US', {
                year: 'numeric',
                month: 'short',
                day: 'numeric',
                hour: 'numeric',
                minute: 'numeric',
                hour12: false,
            });
    }
    if (typeof frontmatter.location_url === 'array' && frontmatter.location_url.length == 1) {
        frontmatter.location_url = frontmatter.location_url[0];
    }

    let event_location = '';
    if (typeof frontmatter.location_url === 'string') {
        event_location = frontmatter.location_url;
    } else {
        event_location = frontmatter.location_url[0];
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

<div class={'card mb-2 rounded-0 rounded-end border-' + type_class}>
    <div class="card-body">
        <div class={'card-title'}>
            <a href={'events/' + slug}>
                {frontmatter.title}
            </a>
            {#if time_category === 'current'}
                <div class="dropwdown btn-group float-end" role="group">
                    <button
                        type="button"
                        class="btn btn-outline-secondary dropdown-toggle"
                        href="#"
                        data-bs-toggle="dropdown"
                        aria-haspopup="true"
                        aria-expanded="false"
                    >
                        Watch now
                    </button>
                    <div class="dropdown-menu text-secondary">
                        {#if typeof frontmatter.location_url === 'string'}
                            <a class="dropdown-item" href={frontmatter.location_url} target="_blank" rel="noreferrer">
                                {frontmatter.location_url}
                            </a>
                        {:else}
                            {#each frontmatter.location_url as url}
                                <a class="dropdown-item" href={url} target="_blank" rel="noreferrer"> {url} </a>
                            {/each}
                        {/if}
                    </div>
                </div>
            {/if}
        </div>
        <div class="card-text">
            <p>{frontmatter.subtitle}</p>
            <div class="d-flex justify-content-between align-items-end">
                <div>
                    <p class="text-muted">
                        {event_date}
                    </p>
                </div>
                <div class="btn-group" role="group" aria-label="Basic example">
                    <a href={'events/' + slug} class="btn btn-outline-success text-nowrap">See details</a>
                    {#if time_category !== 'past'}
                        <div class="dropwdown btn-group" role="group">
                            <button
                                type="button"
                                class="btn btn-outline-success dropdown-toggle"
                                href="#"
                                data-bs-toggle="dropdown"
                                aria-haspopup="true"
                                aria-expanded="false"
                            >
                                <i class="far fa-calendar-plus me-1" aria-hidden="true" /> Export event
                            </button>
                            <div class="dropdown-menu">
                                <button class="dropdown-item" on:click={() => downloadIcal()}> Download iCal Event</button>
                                <a class="dropdown-item" href={googleCalendar} target="_blank"  rel="noreferrer">
                                    Add to Google Calendar</a
                                >
                                <a class="dropdown-item" href={outlookCalendar} target="_blank"  rel="noreferrer">
                                    Add to Microsoft Outlook</a
                                >
                            </div>
                        </div>
                    {/if}
                </div>
            </div>
        </div>
    </div>
</div>

<style lang="scss">
    @import '../styles/_variables';
    .card .card-title a {
        color: $success;
    }
    .card.rounded-0 {
        border-top: 0;
        border-right: 0;
        border-bottom: 0;
        border-left: 5px solid;
        // overflow: hidden;
    }
</style>
