<script lang="ts">
    import { newBlogPost } from '@components/store';
</script>

{#if newBlogPost}
    <span
        class="float-end text-primary-emphasis border border-primary badge rounded-pill fw-normal text btn-sm text-nowrap"
    >
        new post
    </span>
{/if}
