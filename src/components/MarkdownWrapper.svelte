<script lang="ts">
    import { onMount } from 'svelte';
    import { currentHeading } from '@components/store';

    export let headings: { text: string; slug: string; depth: number; fa_icon?: string }[] = [];

    // find current heading in viewport with IntersectionObserver
    onMount(() => {
        const observer = new IntersectionObserver(
            (entries) => {
                entries.forEach((entry) => {
                    if (entry.isIntersecting) {
                        currentHeading.set(entry.target.id);
                    }
                });
            },
            {
                // root: document.querySelector('.markdown-content'),
                rootMargin: '0px 0px -92% 0px',
                // rootMargin: '100px 0px -' + activeHeadingScrollOffset + 'px 0px',
            }
        );
        headings.forEach((heading) => {
            const element = document.querySelector('#' + heading.slug);
            observer.observe(element);
        });
    });
</script>

<div class="markdown-content">
    <slot />
</div>
