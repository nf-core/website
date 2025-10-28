<script lang="ts">
    import { showHelp } from "@components/store";
    import { slide } from "svelte/transition";

    interface Props {
        buttonText?: string;
        buttonClass?: string;
        textClass?: string;
        children?: import("svelte").Snippet;
    }

    let {
        buttonText = "Help text",
        buttonClass = "btn-outline-secondary me-auto",
        textClass = "",
        children,
    }: Props = $props();

    let isOpen = $state(false);
    let userClosed = $state(false); // Track if user explicitly closed it

    // If user closed it, stay closed. Otherwise follow showHelp or isOpen
    let open = $derived(!userClosed && ($showHelp || isOpen));

    function toggleShow() {
        if (open) {
            // User is closing it
            isOpen = false;
            userClosed = true;
        } else {
            // User is opening it
            isOpen = true;
            userClosed = false;
        }
    }
</script>

<div>
    <div class="d-flex">
        <button class={"btn " + buttonClass} class:open type="button" onclick={toggleShow}>
            {buttonText}
        </button>
    </div>
    {#if open}
        <div
            class={"collapse collapsible-content " + textClass}
            class:show={open}
            transition:slide={{ duration: 300, axis: "y" }}
            style="padding-top: 1px;"
        >
            {@render children?.()}
        </div>
    {/if}
</div>

<style lang="scss">
    .btn.open {
        border-bottom-left-radius: 0;
        border-bottom-right-radius: 0;
        transition:
            border-bottom-left-radius 0.3s ease-out,
            border-bottom-right-radius 0.3s ease-out;
        border-bottom: 0;
    }
    .btn:not(.open) {
        transition: border-bottom 0.1s ease 0.25s; // Delay restoring border
    }
</style>
