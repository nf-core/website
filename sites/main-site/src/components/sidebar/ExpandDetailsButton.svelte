<script lang="ts">
    let expanded = $state(false);

    const toggleDetails = () => {
        expanded = !expanded;

        if (expanded) {
            document.querySelectorAll("details:not([open])").forEach((details) => {
                details.setAttribute("open", "");
            });
        } else {
            document.querySelectorAll("details[open]").forEach((details) => {
                details.removeAttribute("open");
            });
        }
    };
</script>

<button class="btn btn-sm btn-outline-secondary expand-details w-100" onclick={toggleDetails}>
    <span>
        {#if expanded}
            <i class="fa-solid fa-arrows-to-line"></i>
            Hide details
        {:else}
            <i class="fa-solid fa-arrows-from-line"></i>
            Show details
        {/if}
    </span>
</button>
