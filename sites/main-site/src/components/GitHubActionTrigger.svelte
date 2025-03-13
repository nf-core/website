<script>
    let isLoading = $state(false);
    let result = $state(null);

    async function triggerGitHubAction() {
        isLoading = true;
        result = null;

        try {
            const response = await fetch("/api/trigger-github-action", {
                method: "POST",
            });

            if (response.ok) {
                result = { success: true, message: "GitHub Action triggered successfully!" };
            } else {
                result = { success: false, message: "Failed to trigger GitHub Action" };
            }
        } catch (error) {
            console.error("Error:", error);
            result = { success: false, message: "Error triggering GitHub Action" };
        } finally {
            isLoading = false;
        }
    }
</script>

<div class="github-action-trigger">
    <button onclick={() => triggerGitHubAction()} disabled={isLoading} class:loading={isLoading}>
        {#if isLoading}
            Triggering...
        {:else}
            Trigger GitHub Action
        {/if}
    </button>

    {#if result}
        <div class="result" class:success={result.success} class:error={!result.success}>
            {result.message}
        </div>
    {/if}
</div>
