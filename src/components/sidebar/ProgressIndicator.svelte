<script>
    import { onMount } from 'svelte';

    export let progress = 0;

    onMount(() => {});
</script>

<div class="progress-bar">
    <div class="progress" style="width: {progress}%"></div>
</div>

<style>
    .progress-bar {
        width: 100%;
        height: 2px;
        background-color: var(--bs-light);
    }

    .progress-bar .progress {
        height: 100%;
        background-color: var(--bs-success);
        transition: width 0.3s;
    }
</style>
