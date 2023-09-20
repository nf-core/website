<script lang="ts">
    import VideoButton from '@components/VideoButton.svelte';
    import ExportEventButton from '@components/event/ExportEventButton.svelte';
    import { onMount } from 'svelte';

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
    export let showDescription: boolean = true;
    export let narrow: boolean = false;

    const event_duration = (event) => {
        let duration;
        if (event.start_date === event.end_date) {
            duration =
                new Date(event.start).toLocaleString('en-US', {
                    year: 'numeric',
                    month: 'short',
                    day: 'numeric',
                    hour: 'numeric',
                    minute: 'numeric',
                    hour12: false,
                }) +
                '-' +
                new Date(event.end).toLocaleString('en-US', {
                    hour: 'numeric',
                    minute: 'numeric',
                    hour12: false,
                });
        } else {
            duration =
                new Date(event.start).toLocaleString('en-US', {
                    year: 'numeric',
                    month: 'short',
                    day: 'numeric',
                    hour: 'numeric',
                    minute: 'numeric',
                    hour12: false,
                }) +
                '<wbr> - <wbr>' +
                new Date(event.end).toLocaleString('en-US', {
                    year: 'numeric',
                    month: 'short',
                    day: 'numeric',
                    hour: 'numeric',
                    minute: 'numeric',
                    hour12: false,
                });
        }
        return duration;
    };
    const event_type_classes = {
        bytesize: 'success',
        hackathon: 'primary',
        talk: 'info',
        training: 'warning',
    };

    let event_date = event_duration(frontmatter);

    const type_class = event_type_classes[type];
    onMount(() => {
        event_date = event_duration(frontmatter);
    });
</script>

<div class={'card mb-3 rounded-0 rounded-end ' + type} style="border-left-color:var(--bs-{type_class});">
    <div class="card-body">
        <div class="card-title">
            <h4 id={'event-' + slug.split('/')[1]}>
                <a class="text-center" href={/events/ + slug + '/'}>
                    {frontmatter.title}
                </a>
                {#if time_category === 'current' && frontmatter.location_url}
                    <div class="float-end d-none d-md-inline">
                        <VideoButton urls={frontmatter.location_url} btnClass="btn-danger" />
                    </div>
                {/if}
            </h4>
        </div>
        <div class="card-text">
            {#if showDescription}
                <p class="mb-0">{frontmatter.subtitle}</p>
            {/if}
            <div
                class="d-flex align-items-center mt-2 flex-wrap justify-content-start"
                class:justify-content-md-end={!narrow}
            >
                <p class="text-nowrap text-center text-md-start pe-3 mt-2 ms-1" class:d-md-none={!narrow}>
                    <i class="fa-regular fa-calendar me-2" />{@html event_date}
                </p>
            </div>
        </div>
    </div>
    <div class="card-footer p-0" class:p-md-2={!narrow}>
        <div class="d-flex align-items-center justify-content-between">
            <p class="d-none text-wrap mb-0 ms-2 align-middle" class:d-md-inline-block={!narrow}>
                {@html event_date}
            </p>
            <div
                class="btn-group float-end"
                class:w-100={narrow}
                class:narrow
                role="group"
                aria-label="See details or export calendar event"
            >
                <a
                    href={'/events/' + slug + '/'}
                    class="btn btn-outline-success text-nowrap rounded-start-0"
                    class:rounded-0={['current', 'future'].includes(time_category)}>See details</a
                >
                {#if time_category === 'future'}
                    <ExportEventButton {frontmatter} add_class={'btn-outline-success ' + ' rounded-top-0'} />
                {/if}
            </div>
        </div>
        {#if time_category === 'current' && frontmatter.location_url}
            <VideoButton
                urls={frontmatter.location_url}
                btnClass=" d-md-none btn-danger w-100 rounded-top-0 rounded-start-0"
            />
        {/if}
    </div>
</div>

<style lang="scss">
    @import '@styles/_variables.scss';
    .card.rounded-0 {
        border-left: 5px solid;
    }
    .narrow .btn:first-child {
        border-left: 0;
    }
    @include media-breakpoint-up(md) {
        .btn-group.float-end:not(.narrow) {
            .btn:first-child {
                border-top-left-radius: $border-radius !important;
                border-bottom-left-radius: $border-radius !important;
            }
            :global(.btn.dropdown-toggle) {
                border-top-right-radius: $border-radius !important;
                border-bottom-right-radius: $border-radius !important;
            }
        }
    }
    @include media-breakpoint-down(md) {
        .btn-group.float-end {
            width: 100%;
            .btn:first-child {
                border-left: 0;
            }
        }
    }
</style>
