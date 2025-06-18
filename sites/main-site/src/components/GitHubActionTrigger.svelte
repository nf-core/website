<script>
    import { Octokit } from "octokit";
    export let pipelineName = ""; // Optional pipeline name to update

    let isLoading = false;
    let result = null;
    let alreadyRunning = false;

    async function triggerPipelineUpdate() {
        isLoading = true;
        result = null;
        alreadyRunning = false;

        const octokit = new Octokit({
            auth: import.meta.env.GITHUB_TOKEN,
        });

        try {
            // First check if workflow is already running
            const { data: runningWorkflows } = await octokit.rest.actions.listWorkflowRuns({
                owner: "nf-core",
                repo: "website",
                workflow_id: "build-json-files.yml",
                status: "in_progress",
                event: "workflow_dispatch",
            });
            console.log(runningWorkflows);

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
            await octokit.rest.repos.createDispatchEvent({
                owner: "nf-core",
                repo: "website",
                event_type: "update-pipeline",
                ref: "main",
                client_payload: {
                    pipeline_name: pipelineName,
                },
            });

            result = { success: true, message: "Pipeline update triggered successfully!" };
        } catch (error) {
            console.error("Error:", error);
            if (error.status === 401) {
                result = {
                    success: false,
                    message: "Authentication failed - please check the GitHub token configuration",
                };
            } else {
                result = {
                    success: false,
                    message: `Error triggering pipeline update: ${error.message}`,
                };
            }
        } finally {
            isLoading = false;
        }
    }
</script>

<div class="pipeline-trigger d-flex align-items-center ms-2">
    <button
        onclick={triggerPipelineUpdate}
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
        <div class="result small" class:success={result.success} class:error={!result.success}>
            {result.message}
            {#if alreadyRunning && result.url}
                <a href={result.url} target="_blank" rel="noopener noreferrer">View running workflow</a>
            {/if}
        </div>
    {/if}
</div>

<style>
</style>
