<script lang="ts">
    export let text: string;
    $: copied = false;
    const copyToClipboard = (text: string) => {
        navigator.clipboard.writeText(text);

        copied = true;

        setTimeout(() => {
            copied = false;
        }, 1000);
        return true;
    };
</script>

<span class="copy-url" on:click={() => copyToClipboard(text)} on:keypress={() => copyToClipboard(text)}>
    <slot />
    {#if copied}
        <i class="fa-regular px-1 fa-clipboard-check" />
    {/if}
</span>
