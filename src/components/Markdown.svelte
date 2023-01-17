<script>
    import { onMount } from 'svelte';
    export let headings = [];
    function set_active_heading() {
        let idx = headings.length;
        while (idx--) {
            const { top } = headings[idx].getBoundingClientRect();
            // loop through headings from last to first until we find one that the viewport already
            // scrolled past. if none is found, set make first heading active
            if (top < activeHeadingScrollOffset || idx === 0) {
                activeHeading = headings[idx];
                activeTocLi = tocItems[idx];
                if (keepActiveTocItemInView && activeTocLi) {
                    // get the currently active ToC list item
                    // scroll the active ToC item into the middle of the ToC container
                    nav.scrollTo?.({ top: activeTocLi?.offsetTop - nav.offsetHeight / 2 });
                }
                return; // exit while loop if updated active heading
            }
        }
    }
</script>
<div class="markdown-content" on:scroll={set_active_heading}>
    <slot />
</div>

<style lang="scss">
    .markdown-content pre.astro-code {
        padding: 0;
    }
</style>
