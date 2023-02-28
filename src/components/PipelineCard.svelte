<script>
    import { formatDistanceToNow } from 'date-fns';
    export let pipeline;
    const name = pipeline.name;
    const body = pipeline.description;
    const stars = pipeline.stargazers_count;
    const topics = pipeline.topics;
    const releases = pipeline.releases;
    const archived = pipeline.archived;
    const released = releases.length > 1;
    var latest_release, tag_name, release_date_ago;
    if (released) {
        latest_release = releases[0];
        tag_name = latest_release.tag_name;
        release_date_ago = formatDistanceToNow(new Date(latest_release.published_at));
    }
    const href = pipeline.name + '/' + (released ? tag_name : 'dev');
</script>

<div class="card flex-fill m-2">
    <div class="card-header border-bottom-0 bg-transparent">
        <h2 class="mb-0 d-flex justify-content-between align-items-center">
            <a {href}
                >{name}
                {#if archived}
                    <i class="fa-solid fa-archive text-info" />
                {:else if released}
                    <i class="fa-solid fa-check text-success" title="released" data-bs-toggle="tooltip" />
                {:else}
                    <i class="fa-solid fa-wrench text-warning" />
                {/if}
            </a>
            <small class="gh-stats text-small">
                <span>
                    {#if released}
                        <a
                            href={'https://github.com/nf-core/' + name + '/releases/tag/' + tag_name}
                            style={{ cursor: 'pointer' }}
                            class="text-muted text-decoration-none"
                        >
                            <i class="fa-regular fa-tag ms-3 me-1" />
                            {tag_name}
                        </a>
                    {/if}
                </span>
                <a
                    href={'https://github.com/nf-core/' + name + '/stargazers'}
                    target="_blank"
                    rel="noreferrer"
                    class="stargazers text-decoration-none mt-2 ms-2 text-warning"
                    title=""
                    data-bs-toggle="tooltip"
                    data-html="true"
                    data-bs-original-title={stars + ' stargazers on GitHub'}
                    style={{ cursor: 'pointer' }}
                >
                    <i class="fa-regular fa-star" aria-hidden="true" />

                    {stars}
                </a>
            </small>
        </h2>
    </div>
    <div class="card-body pt-0 d-flex flex-column">
        <p class="topics mt-0 mb-0">
            {#each topics as topic}
                <span class="badge bg-body-tertiary text-success me-2">{topic}</span>
            {/each}
        </p>
        {#if body}
            <p class="description flex-grow-1 mb-0">{body}</p>
        {/if}

        {#if released}
            <p class="text-muted align">Last release {release_date_ago}</p>
        {/if}
    </div>
</div>

<style>
    p {
        margin-top: 0.5rem;
    }
    /* .link-card:is(:hover, :focus-within) {
        background-position: 0;
    }
    .link-card:is(:hover, :focus-within) h2 {
        color: rgb(var(--accent));
    } */
    .card {
        max-width: 40rem;
    }
    .badge.text-success {
        font-weight: 400;
    }
    .gh-stats {
        float: right;
    }
    .gh-stats a:hover .fa-regular {
        font-weight: 900;
    }
</style>
