<script>
    export let headings = [];
    // get minimal heading depth from headings
    let max_heading_depth = 2
    headings = headings.filter(h => h.depth <= max_heading_depth);

</script>

<!-- Render remark-toc as a nav-list -->
<div class="nav flex-column position-fixed">
    <strong class="d-none d-md-block h6 my-2 text-body">On this page</strong>
    <hr class="d-none d-md-block my-1" />

    <nav id="TableOfContents">
        <ul>
            {#each headings as heading}
            <li class={"nav-item ps-"+heading.depth}>
                <a class="nav-link small py-1" href={"#"+heading.slug}>
                    {heading.text}
                </a>
            </li>
            {/each}
        </ul>
    </nav>
</div>

<style lang="scss">
    nav>ul{
        list-style: none;
        padding-left: 0;
    }
</style>
