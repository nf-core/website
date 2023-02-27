<script lang="ts">
    import { onMount } from 'svelte';

    export let cmd: string;

    let Tooltip: any;
    $: copied = false;

    onMount(() => {
        // activate bootstrap tooltips
        const tooltipTriggerList = document.querySelectorAll('.copy-txt');
        const tooltipElement = [...tooltipTriggerList][0];
        Tooltip = bootstrap.Tooltip.getInstance(tooltipElement);
        if (!Tooltip) {
            Tooltip = new bootstrap.Tooltip(tooltipElement);
        }
    });

    const copyToClipboard = (text: string) => {
        navigator.clipboard.writeText(text);
        copied = true;
        Tooltip.setContent({ '.tooltip-inner': 'Copied!' });
        debugger;
        setTimeout(() => {
            copied = false;

            Tooltip.hide();
            Tooltip.setContent({ '.tooltip-inner': 'Copy to clipboard' });
        }, 1000);
        return true;
    };
</script>

<div class="input-group module-install-cmd">
    <span class="input-group-text bg-body border-secondary border-end-0"><i class="fa-regular fa-terminal" /></span>
    <input
        type="text"
        class="form-control input code bg-body border-secondary border-start-0"
        id="module-install-cmd-text"
        data-autoselect=""
        value={cmd}
        readonly={true}
        on:click={() => copyToClipboard(cmd)}
    />
    <button
        class="btn btn-outline-secondary bg-body-secondary copy-txt"
        data-bs-target="#module-install-cmd-text"
        data-bs-toggle="tooltip"
        data-bs-placement="left"
        title="Copy to clipboard"
        type="button"
        on:click={() => copyToClipboard(cmd)}
        ><i class={'fa-regular px-1 fa-clipboard' + (copied ? '-check' : '')} /></button
    >
</div>

<style lang="scss">
    @import '../styles/_variables.scss';

    input:focus {
        box-shadow: none;
    }

    @include color-mode(dark) {
        .border-secondary {
            border-color: $gray-900 !important;
        }
    }
    .copy-txt {
        // border-width: pt;
        border-color: $secondary;
    }
</style>
