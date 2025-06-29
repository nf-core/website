<script lang="ts">
    import ListingCard from "@components/ListingCard.svelte";
    import TagSection from "@components/TagSection.svelte";
    import Markdown from "@components/markdown/Markdown.svelte";
    import { formatDistanceToNow, add } from "date-fns";
    import { Confetti } from "svelte-confetti";

    let { pipeline } = $props();

    const name = pipeline.name;
    const body = pipeline.description;
    const stars = pipeline.stargazers_count;
    const topics = pipeline.topics;
    const releases = pipeline.releases;
    const archived = pipeline.archived;
    const released = releases.length > 1;
    let latestRelease,
        tagName = $state(),
        releaseDateAgo = $state(),
        recentRelease = $state(false);
    const lastChangesAge = formatDistanceToNow(new Date(releases[0].published_at), {
        addSuffix: true,
    });
    if (released) {
        latestRelease = releases[0];
        tagName = latestRelease.tag_name;
        releaseDateAgo = formatDistanceToNow(new Date(latestRelease.published_at), {
            addSuffix: true,
        });
        // Check if release is less than 1 day old
        recentRelease =
            new Date(latestRelease.published_at).getTime() > add(new Date().getTime(), { days: -1 }).getTime();
    }
</script>

<ListingCard {recentRelease}>
    {#snippet cardHeader()}
        <div class="d-flex flex-wrap justify-content-between align-items-center">
            <a class="text-decoration-none" href={"/" + pipeline.name + "/" + (released ? tagName : "dev") + "/"}
                >{name}
                {#if archived}
                    <i class="fa-solid fa-xs fa-archive text-info" title="archived" data-bs-toggle="tooltip"></i>
                {:else if released}
                    <i class="fa-solid fa-xs fa-check text-success" title="released" data-bs-toggle="tooltip"></i>
                {:else}
                    <i class="fa-solid fa-xs fa-wrench text-warning" title="under development" data-bs-toggle="tooltip"
                    ></i>
                {/if}
            </a>

            <small class="gh-stats fs-5">
                <a
                    href={"https://github.com/nf-core/" + name + "/stargazers"}
                    target="_blank"
                    rel="noreferrer"
                    class="stargazers text-decoration-none mt-2 ms-2 text-warning"
                    title={stars + " stargazers on GitHub"}
                    data-bs-toggle="tooltip"
                    data-html="true"
                    data-bs-original-title={stars + " stargazers on GitHub"}
                >
                    <i class="fa-regular fa-star" aria-hidden="true"></i>

                    {stars}
                </a>
            </small>
        </div>
    {/snippet}

    {#snippet cardBody()}
        <div class="d-flex flex-column justify-content-between h-100">
            <div class="recent-release-badge text-center">
                {#if recentRelease}
                    <a
                        href={"https://github.com/nf-core/" + name + "/releases/tag/" + tagName}
                        class="text-decoration-none badge text-bg-success fs-6 mb-0 rounded-top-0"
                        >New release! <Confetti x={[-1.5, 1.75]} amount={100} rounded={true} /></a
                    >
                {/if}
            </div>
            {#if body}
                <div class="description flex-grow-1" class:pt-1={recentRelease}><Markdown md={body} /></div>
            {/if}

            <div class="mb-2">
                <TagSection tags={topics} type="keywords" />
            </div>

            {#if released}
                <p class="release">
                    <a
                        role="button"
                        href={"https://github.com/nf-core/" + name + "/releases/tag/" + tagName}
                        class="btn btn-outline-secondary"
                    >
                        <i class="fa-regular fa-tag me-1"></i>
                        {tagName}
                    </a>
                    <span class="text-body-secondary text-small"> released {releaseDateAgo}</span>
                </p>
            {:else}
                <p class="release">
                    <span class="text-body-secondary text-small"> last changes {lastChangesAge}</span>
                </p>
            {/if}
        </div>
    {/snippet}
</ListingCard>

<style lang="scss">
    .gh-stats {
        float: right;
    }
    .gh-stats a,
    .release a {
        &:hover {
            text-decoration: underline !important;
            .fa-regular {
                font-weight: 900;
            }
        }
    }
    .recent-release-badge {
        margin-top: -1px;
    }
</style>
