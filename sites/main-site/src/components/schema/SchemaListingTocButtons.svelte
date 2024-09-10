<script lang="ts">
    import { showHidden, showHelp } from '@components/store';
</script>

<div class="btn-group btn-sm w-100 px-3 px-xl-0" role="group" aria-label="Show/Hide hidden parameters and help texts">
    <button type="button" class="btn btn-outline-secondary btn-sm" on:click={() => showHidden.set(!$showHidden)}
        >{$showHidden ? 'Hide' : 'Show'} hidden</button
    >
    <button type="button" class="btn btn-outline-secondary btn-sm" on:click={() => showHelp.set(!$showHelp)}
        >{$showHelp ? 'Hide' : 'Show'} all help texts</button
    >
</div>
