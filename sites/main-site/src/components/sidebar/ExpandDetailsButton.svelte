<script lang="ts">
    let expanded = false;

    function ToggleDetails() {
        expanded = !expanded;
        if (expanded) {
            document.querySelectorAll('details:not([open])').forEach((details) => {
                details.setAttribute('open', '');
            });
        } else {
            document.querySelectorAll('details[open]').forEach((details) => {
                details.removeAttribute('open');
            });
        }
    }
</script>

<button class="btn btn-sm btn-outline-secondary expand-details w-100" on:click={ToggleDetails}>
    <span class="">
        {#if expanded}
            <i class="fa-solid fa-arrows-to-line" />
            Hide details
        {:else}
            <i class="fa-solid fa-arrows-from-line" />
            Show details
        {/if}
    </span></button
>
