<script lang="ts">
    export let timeSpan: [number, number] = [0, 0];

    const now = new Date().getTime();
    const show = timeSpan[0] < now && timeSpan[1] > now;
</script>

{#if show}
    <span
        class="float-end text-primary-emphasis border border-primary badge rounded-pill fw-normal text btn-sm text-nowrap"
    >
        new post
    </span>
{/if}
