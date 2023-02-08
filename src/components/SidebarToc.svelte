<script lang="ts">
    export let headings: {
        text: string;
        slug: string;
        depth: number;
        fa_icon?: string;
    }[];

    // get minimal heading depth from headings
    let max_heading_depth = 3;
    // filter out headings that are higher than max_heading_depth
    headings = headings.filter((h) => h.depth <= max_heading_depth);

    let min_heading_depth = Math.min(...headings.map((h) => h.depth));
    // make margin classes from min to max heading depth
    let headingMargin = {};
    for (let i = min_heading_depth; i <= 6; i++) {
        headingMargin[i] = 'ps-' + (i - min_heading_depth) * 2;
    }
</script>

<div class="nav flex-column sticky-top-under">
    <div class="d-none d-md-block">
        <strong class="h6 my-2 text-body">On this page</strong>
        <hr class="my-1" />
        <nav id="TableOfContents d-none d-md-flex flex-column ">
            <ul class="mb-0">
                {#each headings as heading (heading)}
                    <li class={'nav-item ' + headingMargin[heading.depth]} class:active={heading.depth==2}>
                        <a class="nav-link py-1 ps-1 text-muted" href={'#' + heading.slug}>
                            {#if heading.fa_icon}
                                <i class={heading.fa_icon} aria-hidden="true" />
                            {/if}
                            {@html heading.text}
                        </a>
                    </li>
                {/each}
            </ul>
            <div class="">
                <a
                    href="#/"
                    class="back-to-top text-muted text-small float-end mb-2"
                    on:click={() => window.scrollTo(0, 0)}
                >
                    <i class="fa-solid fa-arrow-up-to-line" aria-hidden="true" /> Back to top
                </a>

                <slot />
            </div>
        </nav>
    </div>
</div>

<style lang="scss">
    @import 'src/styles/_variables.scss';
    nav > ul {
        font-size: 0.875rem;
        list-style: none;
        padding-left: 0;
        overflow-y: auto;
        max-height: calc(100% - 56rem);
    }
    .sticky-top-under{
        top: 5.05rem;
        position: sticky;
    }
    li{
        border-left: 2pt solid $secondary;
    }
    li.active{
        border-left: 2pt solid $success;
        background-color: transparentize($success, 0.75);
    }
</style>
