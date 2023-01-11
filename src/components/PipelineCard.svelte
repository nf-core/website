<script>
    import { format } from 'timeago.js';
    export let pipeline;
    const href = pipeline.name;
    const name = pipeline.name;
    const body = pipeline.description;
    const stars = pipeline.stargazers_count;
    const topics = pipeline.topics;
    const releases = pipeline.releases;
    const archived = pipeline.archived;
    const released = releases.length > 0;
    var latest_release, tag_name, release_date_ago;

    if (released) {
        latest_release = releases[releases.length - 1];
        tag_name = latest_release.tag_name;
        release_date_ago = format(new Date(latest_release.published_at), 'en_GB');
    }
</script>

<div class="card w-100 p-3 pb-2 m-2">
    <div class="card-name ">
        <h2 class="mb-0 d-flex justify-content-between align-items-center">
            <a {href}>{name}
            {#if archived}
                <i class="fa-solid fa-archive text-info" />
            {:else if released}
                <i class="fa-solid fa-check text-success" title="released" data-bs-toggle="tooltip"/>
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
    <div class="card-body py-0 d-flex flex-column">
        <p class="topics mt-0 mb-0">
            {#each topics as topic}
                <span class="badge bg-body-tertiary text-success mx-1">{topic}</span>
            {/each}
        </p>
        <p class="description flex-grow-1 mb-0">{body}</p>

    </div>
    {#if released}
            <p class="text-muted align">Last release {release_date_ago}</p>
        {/if}
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
        max-width: 40%;
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
