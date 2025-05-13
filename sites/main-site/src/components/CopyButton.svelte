<script lang="ts">
    interface Props {
        text: string;
        label?: string;
        copiedLabel?: string;
        classes?: string;
        copiedClasses?: string;
        children?: import("svelte").Snippet;
    }

    let {
        text,
        label = "",
        copiedLabel = '<i class="fa-regular px-1 fa-clipboard-check" />',
        classes = "",
        copiedClasses = "",
        children,
    }: Props = $props();

    let copied = $state(false);

    let currentClasses = $derived(copied ? copiedClasses : classes);
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
    class={"copy-url " + currentClasses}
    onclick={() => copyToClipboard(text)}
    onkeypress={() => copyToClipboard(text)}
    data-bs-toggle="tooltip"
    title="Copy to clipboard"
    role="button"
    tabindex="0"
    >{@render children?.()}{#if copied}{@html copiedLabel}
    {:else}{@html label}{/if}</span
>

<style>
    .copy-url {
        /* Set the height of the button,
        because the height is different
        with and without text */
        height: 2rem;
    }
</style>
