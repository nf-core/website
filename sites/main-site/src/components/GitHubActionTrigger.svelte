<script>
    export let pipelineName = ""; // Optional pipeline name to update

    let isLoading = false;
    let result = null;
    let alreadyRunning = false;

    async function triggerPipelineUpdate() {
        isLoading = true;
        result = null;
        alreadyRunning = false;

        try {
            // First check if workflow is already running
            const checkResponse = await fetch(
                "https://api.github.com/repos/nf-core/website/actions/runs?status=in_progress&event=workflow_dispatch",
                {
                    headers: {
                        Accept: "application/vnd.github.v3+json",
                        Authorization: `token ${import.meta.env.GITHUB_TOKEN}`, // Use public env var
                    },
                },
            );

            const runningWorkflows = await checkResponse.json();

            // Check if our specific workflow is already running with the same pipeline
            const duplicateRun = runningWorkflows.workflow_runs?.find(
                (run) =>
                    run.name.toLowerCase().replace(/ /g, "_") === "build_json_files" &&
                    ((!pipelineName && !run.display_title.includes("pipeline_name")) ||
                        (pipelineName && run.display_title.includes(pipelineName))),
            );

            if (duplicateRun) {
                alreadyRunning = true;
                result = {
                    success: false,
                    message: `This workflow is already running. Started at ${new Date(duplicateRun.created_at).toLocaleTimeString()}`,
                    url: duplicateRun.html_url,
                };
                return;
            }

            // If not already running, trigger the workflow
            const response = await fetch("https://api.github.com/repos/nf-core/website/dispatches", {
                method: "POST",
                headers: {
                    Accept: "application/vnd.github.v3+json",
                    Authorization: `token ${import.meta.env.PUBLIC_GITHUB_TOKEN}`,
                    "Content-Type": "application/json",
                },
                body: JSON.stringify({
                    event_type: "update-pipeline",
                    client_payload: {
                        pipeline_name: pipelineName,
                    },
                }),
            });

            if (response.status === 204) {
                result = { success: true, message: "Pipeline update triggered successfully!" };
            } else {
                result = { success: false, message: `Failed to trigger pipeline update. Status: ${response.status}` };
            }
        } catch (error) {
            console.error("Error:", error);
            result = { success: false, message: "Error triggering pipeline update" };
        } finally {
            isLoading = false;
        }
    }
</script>

<div class="pipeline-trigger d-flex align-items-center">
    <button
        on:click={triggerPipelineUpdate}
        disabled={isLoading}
        class:loading={isLoading}
        class="btn btn-sm btn-outline-secondary"
    >
        {#if isLoading}
            <i class="fa-solid fa-fw fa-rotate-right fa-spin"></i>
        {:else if pipelineName}
            <i class="fa-solid fa-fw fa-rotate-right"></i>
        {:else}
            Update all pipelines
        {/if}
    </button>

    {#if result}
        <div class="result" class:success={result.success} class:error={!result.success}>
            {result.message}
            {#if alreadyRunning && result.url}
                <a href={result.url} target="_blank" rel="noopener noreferrer">View running workflow</a>
            {/if}
        </div>
    {/if}
</div>

<style>
</style>
