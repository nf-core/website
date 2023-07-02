<script lang="ts">
    import VideoButton from '@components/VideoButton.svelte';
    import ExportEventButton from '@components/event/ExportEventButton.svelte';

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
    export let slug: string = '';
    export let type: string = '';
    export let time_category: string = '';
    export let event_type_classes: {}[] = [];

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
            '<wbr> - <wbr>' +
            frontmatter.end.toLocaleString('en-US', {
                year: 'numeric',
                month: 'short',
                day: 'numeric',
                hour: 'numeric',
                minute: 'numeric',
                hour12: false,
            });
    }
    const type_class = event_type_classes[type];
</script>

<div class={'card mb-3 rounded-0 rounded-end ' + type} style="border-left-color:var(--bs-{type_class});">
    <div class="card-body">
        <div class="card-title">
            <h3>
                <a class="text-center" href={slug + '/'}>
                    {frontmatter.title}
                </a>
                {#if time_category === 'current'}
                    <div class="float-end d-none d-md-inline">
                        <VideoButton urls={frontmatter.location_url} btnClass="btn-danger" />
                    </div>
                {/if}
            </h3>
        </div>
        <div class="card-text">
            <p class="mb-0">{frontmatter.subtitle}</p>
            <div class="d-flex align-items-center mt-2 flex-wrap justify-content-center justify-content-md-end">
                <p class="text-nowrap d-md-none text-center text-md-start pe-3">
                    <i class="fa-regular fa-calendar me-2" />{@html event_date}
                </p>
            </div>
        </div>
    </div>
    <div class="card-footer p-0 p-md-2">
        <div class="d-flex align-items-center justify-content-between">
            <p class="d-none d-md-inline-block text-wrap mb-0 align-middle">
                {@html event_date}
            </p>
            <div class="btn-group float-end" role="group" aria-label="See details or export calendar event">
                <a href={slug + '/'} class="btn btn-outline-success text-nowrap rounded-top-0 rounded-start-0"
                    >See details</a
                >
                {#if time_category === 'future'}
                    <ExportEventButton {frontmatter} add_class={'btn-outline-success ' + ' rounded-top-0'} />
                {/if}
            </div>
        </div>
        {#if time_category === 'current'}
            <VideoButton
                urls={frontmatter.location_url}
                btnClass=" d-md-none btn-danger w-100 rounded-top-0 rounded-start-0"
            />
        {/if}
    </div>
</div>

<style lang="scss">
    @import '@styles/_variables';
    .card .card-title a {
        // color: $success;
    }
    .card.rounded-0:not(.bytesize) {
        border-left: 5px solid;
    }
    @include media-breakpoint-down(md) {
        .btn-group.float-end {
            width: 100%;
        }
    }
</style>
