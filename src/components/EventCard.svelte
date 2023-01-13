<script>
    import ExportEventButton from "./ExportEventButton.svelte";

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

</script>

<div class={'card mb-3 rounded-0 rounded-end'} style="border-left-color:var(--bs-{type_class});">
    <div class="card-body">
        <div class={'card-title'}>
            <a href={slug}>
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
            <div class="d-flex justify-content-center justify-content-md-between align-items-center flex-wrap flex-md-nowrap">
                <div>
                    <p class="text-muted text-nowrap">
                        {@html event_date}
                    </p>
                </div>
                <div class="btn-group ms-1" role="group" aria-label="See details or export calendar event">
                    <a href={slug} class="btn btn-outline-success text-nowrap">See details</a>
                    {#if time_category !== 'past'}
                        <ExportEventButton {frontmatter}/>
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
        border-left: 5px solid;
    }
</style>
