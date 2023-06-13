<script>
    import { onMount } from 'svelte';
    import { flip } from 'svelte/animate';
    import { fade } from 'svelte/transition';
    import contributors_yml from '../contributors.yaml';

    export let contributors = contributors_yml.contributors
        .filter((contributor) => contributor.image_fn)
        .sort((a, b) => 0.5 - Math.random())
        .slice(0, 10);

    // Homepage contributor images fading in and out
    // TODO: No animation yet
    onMount(() => {
        function shuffle_contributors() {
            let removed_contrib = contributors.shift();
            contributors = contributors; // trigger update
            setTimeout(function () {
                contributors.push(removed_contrib);
                shuffle_contributors();
            }, 3000);
        }
        shuffle_contributors();
    });
</script>

<div id="community" class="homepage-usedby">
    <div class="container py-5">
        <h2>
            <a class="btn btn-success float-end d-none d-md-inline" href="/community#organisations">
                See a complete list &raquo;
            </a>
            <a href="/community#organisations">Used by groups all over the world</a>
        </h2>
        <p>The nf-core community is spread all over the globe and includes a large number of contributing users.</p>
        <p>
            <a class="btn btn-success d-inline d-md-none" href="/community#organisations">
                See a complete list &raquo;
            </a>
        </p>

        <div class="homepage_contrib_logos">
            {#each contributors as contributor (contributor)}
                <a
                    href="/community#{contributor.full_name.toLowerCase().replace(/[^a-z]+/i, '-')}"
                    transition:fade={{ duration: 500 }}
                    animate:flip={{ duration: 500 }}
                >
                    <img
                        src="/src/assets/contributors/white/{contributor.image_fn}"
                        class="my-2 my-lg-3 mx-2"
                        data-bs-placement="bottom"
                        data-bs-toggle="tooltip"
                        title={contributor.full_name}
                        alt={contributor.full_name}
                    />
                </a>
            {/each}
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
            content: '';
            background: #000;
            position: absolute;
            top: 0;
            left: 0;
            bottom: 0;
            right: 0;
            z-index: -2;
        }
        &::after {
            content: '';
            background: url('/flowcell.jpg') no-repeat;
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
        @media (max-width: 576px) {
            img {
                height: 50px;
                margin: 0.5rem 1rem 0.5rem 0;
            }
        }
    }
    .homepage_contrib_logos {
        overflow: hidden;
        white-space: nowrap;
        a {
            display: inline-block;
        }
    }
</style>
