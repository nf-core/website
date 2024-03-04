<script lang="ts">
    export let cmd: string;
    export let icon: boolean = true;
    export let flatTop: boolean = false;
    export let flatBottom: boolean = false;

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

<div class="input-group module-install-cmd">
    {#if icon}
        <span class="input-group-text bg-body border-end-0"><i class="fa-regular fa-terminal" /></span>
    {/if}
    <input
        type="text"
        class="form-control input code bg-body overflow-x-auto"
        class:border-start-0={icon}
        class:rounded-top-0={flatTop}
        class:rounded-bottom-0={flatBottom}
        id="module-install-cmd-text"
        data-autoselect=""
        value={cmd}
        readonly={true}
        on:click={() => copyToClipboard(cmd)}
    />
    <button
        class="btn btn-secondary border copy-txt"
        class:text-bg-success={copied}
        class:rounded-top-0={flatTop}
        data-bs-toggle="tooltip"
        data-bs-placement="left"
        title="Copy to clipboard"
        type="button"
        on:click={() => copyToClipboard(cmd)}
        ><i class={'fa-regular px-1 fa-clipboard' + (copied ? '-check' : '')} /></button
    >
</div>

<style lang="scss">
    @import '@styles/_variables.scss';
    input {
        border-color: $border-color;
        &:focus {
            box-shadow: none;
        }
    }
    .module-install-cmd {
        margin-bottom: -1px;
    }
    .btn.copy-txt {
        background-color: $border-color;
        color: $body-color;
        border-top-right-radius: 3px;
        &:hover {
            background-color: $secondary;
            color: $white;
        }
        :global(.sidebar) & {
            border-bottom-right-radius: 0;
        }
    }

    :global([data-bs-theme='dark']) {
        input {
            border-color: $border-color-dark;
        }
        & .btn.copy-txt {
            background-color: $border-color-dark;
            color: $body-color-dark;
            &:hover,
            &.active {
                background-color: $secondary;
                color: $white;
            }
        }
    }
</style>
