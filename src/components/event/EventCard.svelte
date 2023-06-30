<script lang="ts">
    import ExportEventButton from '@components/event/ExportEventButton.svelte';
    import VideoButton from '@components/VideoButton.svelte';

    export let frontmatter = {
        title: '',
        subtitle: '',
        start: new Date(),
        start_date: new Date(),
        end: new Date(),
        end_date: new Date(),
        type: '',
        location_url: [''],
    };
    export let slug = '';
    export let type = '';
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
            ' - <wbr>' +
            frontmatter.end.toLocaleString('en-US', {
                year: 'numeric',
                month: 'short',
                day: 'numeric',
                hour: 'numeric',
                minute: 'numeric',
                hour12: false,
            });
    }
    const event_type_classes = {
        bytesize: 'success',
        hackathon: 'primary',
        poster: 'danger',
        talk: 'success',
        tutorial: 'info',
        training: 'warning',
    };
    const type_class = event_type_classes[type];
</script>

<div class={'card mb-3 rounded-0 rounded-end'} style="border-left-color:var(--bs-{type_class});">
    <div class="card-body">
        <div class={'card-title'}>
            <a href={slug + '/'}>
                {frontmatter.title}
            </a>
            {#if time_category === 'current'}
                <div class="float-end">
                    <VideoButton urls={frontmatter.location_url} btnClass="btn-danger" />
                </div>
            {/if}
        </div>
        <div class="card-text">
            <p class="mb-2 mb-md-3">{frontmatter.subtitle}</p>
            <div
                class="d-flex justify-content-center justify-content-md-between align-items-center flex-wrap flex-md-nowrap"
            >
                <div>
                    <p class="text-body-secondary text-nowrap text-center text-md-start">
                        {@html event_date}
                    </p>
                </div>
                <div class="btn-group ms-1" role="group" aria-label="See details or export calendar event">
                    <a href={slug + '/'} class="btn btn-outline-success text-nowrap">See details</a>
                    {#if time_category !== 'past'}
                        <ExportEventButton {frontmatter} />
                    {/if}
                </div>
            </div>
        </div>
    </div>
</div>

<style lang="scss">
    @import '@styles/_variables';
    .card .card-title a {
        color: $success;
    }
    .card.rounded-0 {
        border-left: 5px solid;
    }
</style>
