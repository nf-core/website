<script lang="ts">
    export let text: string;
    export let label: string = '';
    export let copiedLabel: string = '<i class="fa-regular px-1 fa-clipboard-check" />';
    export let classes: string = '';
    export let copiedClasses: string = '';

    $: copied = false;
    $: currentClasses = copied ? copiedClasses : classes;
    const copyToClipboard = (text: string) => {
        navigator.clipboard.writeText(text);
        copied = true;
        setTimeout(() => {
            copied = false;
        }, 1000);
        return true;
    };
</script>

<span
    class={'copy-url ' + currentClasses}
    on:click={() => copyToClipboard(text)}
    on:keypress={() => copyToClipboard(text)}
    data-bs-toggle="tooltip"
    title="Copy to clipboard"
    role="button"
    tabindex="0"
    ><slot />{#if copied}
        {@html copiedLabel}
    {:else}
        {@html label}
    {/if}</span
>

<style lang="scss">
    @import '@styles/_variables.scss';
    .copy-url {
        background: $body-bg;
    }
    :global([data-bs-theme='dark']) {
        .copy-url {
            background: $body-bg-dark;
        }
    }
</style>
