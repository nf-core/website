<script>
    import { onMount } from "svelte";
    import contributors_yml from "../config/contributors.yaml";

    const marqueeSpeed = 45;

    let contributors = contributors_yml.contributors
        .filter((contributor) => contributor.image_fn)
        .sort((a, b) => 0.5 - Math.random());

    let displayContributors = $state(contributors.slice(0, 10)); // start by displaying 10 contributors

    function addMoreContributors() {
        let moreContributors = contributors.slice(displayContributors.length, displayContributors.length + 5);
        displayContributors = displayContributors.concat(moreContributors);
    }

    onMount(() => {
        const id = setInterval(addMoreContributors, 10000);
        return () => clearInterval(id);
    });

    const marqueeDurationSec = Math.max(12, Math.round(2000 / marqueeSpeed));
</script>

{#snippet contribStrip(interactive)}
    {#each displayContributors as contributor (contributor)}
        <a
            href="/contributors/#{contributor.full_name.toLowerCase().replace(/[^a-z]+/i, '-')}"
            tabindex={interactive ? undefined : -1}
        >
            <img
                src="/images/contributors/white/{contributor.image_fn}"
                class="my-lg-3 px-2"
                data-bs-placement="bottom"
                data-bs-toggle={interactive ? "tooltip" : undefined}
                title={interactive ? contributor.full_name : undefined}
                alt={interactive ? contributor.full_name : ""}
            />
        </a>
    {/each}
{/snippet}

<div id="community" class="homepage-usedby">
    <div class="container py-5">
        <h2>
            <a class="btn btn-success float-end d-none d-md-inline" href="/contributors#organisations">
                See a complete list &raquo;
            </a>
            <a href="/contributors#organisations">Used by groups all over the world</a>
        </h2>
        <p>The nf-core community is spread all over the globe and includes a large number of contributing users.</p>
        <p>
            <a class="btn btn-success d-inline d-md-none" href="/contributors#organisations">
                See a complete list &raquo;
            </a>
        </p>

        <div class="homepage_contrib_logos">
            <div class="contrib-marquee" style="--marquee-duration: {marqueeDurationSec}s;">
                <div class="contrib-marquee-track">
                    <div class="contrib-marquee-row">
                        {@render contribStrip(true)}
                    </div>
                    <div class="contrib-marquee-row" aria-hidden="true">
                        {@render contribStrip(false)}
                    </div>
                </div>
            </div>
        </div>
    </div>
</div>

<style lang="scss">
    .homepage-usedby {
        position: relative;
        display: block;
        color: #ffffff;
        z-index: 0;
        h2,
        a {
            color: #ffffff;
            text-decoration: none;
        }
        p,
        p a {
            color: rgba(255, 255, 255, 0.8);
        }
        &::before {
            content: "";
            background: #000;
            position: absolute;
            top: 0;
            left: 0;
            bottom: 0;
            right: 0;
            z-index: -2;
        }
        &::after {
            content: "";
            background: url("/flowcell.jpg") no-repeat;
            background-attachment: fixed;
            opacity: 0.8;
            background-size: cover;
            background-position: right top;
            position: absolute;
            top: 0;
            left: 0;
            bottom: 0;
            right: 0;
            z-index: -1;
        }
        img {
            height: 80px;
        }
    }

    .contrib-marquee {
        overflow: hidden;
        width: 100%;
        mask-image: linear-gradient(to right, transparent, #000 4%, #000 96%, transparent);
    }

    .contrib-marquee-track {
        display: flex;
        width: max-content;
        animation: contrib-marquee-scroll var(--marquee-duration, 40s) linear infinite;
    }

    .contrib-marquee:hover .contrib-marquee-track {
        animation-play-state: paused;
    }

    .contrib-marquee-row {
        display: flex;
        flex-shrink: 0;
        align-items: center;
    }

    @keyframes contrib-marquee-scroll {
        from {
            transform: translateX(0);
        }
        to {
            transform: translateX(-50%);
        }
    }
</style>
