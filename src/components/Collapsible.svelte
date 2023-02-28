<script lang="ts">
    import { showHelp } from './store.js';
    export let buttonText: string = 'Help text';
    export let buttonClass: string = 'btn-outline-secondary ms-auto';
    export let textClass: string = '';
    $: show = $showHelp;
</script>

<div>
    <div class="d-flex">
        <button class={'btn ' + buttonClass} class:open={show} type="button" on:click={() => (show = !show)}>
            {buttonText}
        </button>
    </div>
    <div class={'collapse ' + textClass} class:show>
        <slot />
    </div>
</div>

<style lang="scss">
    .btn.open {
        border-bottom-left-radius: 0;
        border-bottom-right-radius: 0;
    }
</style>
