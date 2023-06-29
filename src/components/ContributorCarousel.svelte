<script>
    import Marquee from 'svelte-fast-marquee';
    import contributors_yml from '../config/contributors.yaml';

    let contributors = contributors_yml.contributors
        .filter((contributor) => contributor.image_fn)
        .sort((a, b) => 0.5 - Math.random());
</script>

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
            <Marquee pauseOnHover={true} speed={5}>
                {#each contributors as contributor (contributor)}
                    <a href="/contributors/#{contributor.full_name.toLowerCase().replace(/[^a-z]+/i, '-')}">
                        <img
                            src="/images/contributors/white/{contributor.image_fn}"
                            class="my-lg-3 px-2"
                            data-bs-placement="bottom"
                            data-bs-toggle="tooltip"
                            title={contributor.full_name}
                            alt={contributor.full_name}
                        />
                    </a>
                {/each}
            </Marquee>
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
    }
    :global(.marqueeck-wrapper) {
        background-color: transparent !important;
    }
</style>
