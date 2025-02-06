<script lang="ts">
    import { run } from 'svelte/legacy';

    import { showHelp } from '@components/store';
    interface Props {
        buttonText?: string;
        buttonClass?: string;
        textClass?: string;
        children?: import('svelte').Snippet;
    }

    let {
        buttonText = 'Help text',
        buttonClass = 'btn-outline-secondary me-auto',
        textClass = '',
        children,
    }: Props = $props();
    let show;
    run(() => {
        show = $showHelp;
    });
</script>

<div>
    <div class="d-flex">
        <button class={'btn ' + buttonClass} class:open={show} type="button" onclick={() => (show = !show)}>
            {buttonText}
        </button>
    </div>
    <div class={'collapse ' + textClass} class:show>
        {@render children?.()}
    </div>
</div>

<style lang="scss">
    .btn.open {
        border-bottom-left-radius: 0;
        border-bottom-right-radius: 0;
    }
</style>
