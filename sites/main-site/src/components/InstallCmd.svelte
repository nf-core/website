<script lang="ts">
    interface Props {
        cmd: string;
        icon?: boolean;
        flatTop?: boolean;
        flatBottom?: boolean;
    }

    let { cmd, icon = true, flatTop = false, flatBottom = false }: Props = $props();

    let copied = $state(false);

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
        <span class="input-group-text bg-body border-end-0"><i class="fa-regular fa-terminal"></i></span>
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
        onclick={() => copyToClipboard(cmd)}
    />
    <button
        class="btn btn-secondary border copy-txt"
        class:text-bg-success={copied}
        class:rounded-top-0={flatTop}
        aria-label="Copy to clipboard"
        data-bs-toggle="tooltip"
        data-bs-placement="left"
        title="Copy to clipboard"
        type="button"
        onclick={() => copyToClipboard(cmd)}
        ><i class={"fa-regular px-1 fa-clipboard" + (copied ? "-check" : "")}></i></button
    >
</div>

<style lang="scss">
    input,
    .input-group-text {
        border-color: var(--bs-border-color);
        &:focus {
            box-shadow: none;
        }
    }
    .module-install-cmd {
        margin-bottom: -1px;
    }
    .btn.copy-txt {
        background-color: var(--bs-border-color);
        color: var(--bs-body-color);
        border-top-right-radius: 3px;
        &:hover {
            background-color: var(--bs-secondary);
            color: var(--bs-white);
        }
        :global(.sidebar) & {
            border-bottom-right-radius: 0;
        }
    }
</style>
