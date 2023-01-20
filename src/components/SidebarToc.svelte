<script>
    export let headings = [];
    // get minimal heading depth from headings
    let max_heading_depth = 3;
    // filter out headings that are higher than max_heading_depth
    headings = headings.filter((h) => h.depth <= max_heading_depth);

    let min_heading_depth = Math.min(...headings.map((h) => h.depth));
    // make margin classes from min to max heading depth
    let headingMargin = {};
    for (let i = min_heading_depth; i <= 6; i++) {
        headingMargin[i] = 'ms-' + (i - min_heading_depth) * 2;
    }
</script>

<div class="nav flex-column sticky-top">
    <div class="d-none d-md-block">
        <strong class="h6 my-2 text-body">On this page</strong>
        <hr class="my-1" />
        <nav id="TableOfContents d-none d-md-flex flex-column ">
            <ul>
                {#each headings as heading (heading)}
                    <li class={'nav-item ' + headingMargin[heading.depth]}>
                        <a class="nav-link py-1 ps-1 text-muted" href={'#' + heading.slug}>
                            {heading.text}
                        </a>
                    </li>
                {/each}
            </ul>
            <div class="text-center">
                <button
                    class="btn btn-sm btn-outline-secondary mx-auto back-to-top"
                    on:click={() => window.scrollTo(0, 0)}
                >
                    <i class="fa-solid fa-arrow-up-to-line" aria-hidden="true" /> Back to top
                </button>
            </div>
        </nav>
    </div>
    <!-- dropdown on smaller screens -->

</div>

<style lang="scss">
    @import 'src/styles/_variables';
    .nav {
        padding-top: 5.05rem; // account for navbar
    }
    @include media-breakpoint-down(md) {
        .nav {
            padding-top: 0;
        }
        .toc-md {
            padding-top: 2rem;
            z-index: 9999;
        }
}
    nav > ul {
        font-size: 0.875rem;
        list-style: none;
        padding-left: 0;
        overflow-y: auto;
        max-height: calc(100% - 56rem);
    }
</style>
