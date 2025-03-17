<script lang="ts">
    import ListingCard from "@components/ListingCard.svelte";

    interface Props {
        component: {
            name: string;
            path: string;
            type: string;
            meta: {
                description: string;
                name: string;
                keywords?: string[];
                modules?: string[];
                components?: string[];
            };
            pipelines?: {
                name: string;
                version: string;
            }[];
            subworkflows?: string[];
        };
        pipelines?: {
            name: string;
            version: string;
        }[];
        subworkflows?: string[];
    }

    let { component }: Props = $props();
    const href = "/" + component.type + "s/" + component.name + "/";
    let collapsePipelines = $state(true);
</script>

<ListingCard>
    {#snippet cardHeader()}
        <div>
            <a class="text-decoration-none" {href}>{@html component.name.replace("_", "_<wbr>")} </a>
            <small class="gh-stats text-small"></small>
        </div>
    {/snippet}

    {#snippet cardBody()}
        <div class="d-flex flex-column justify-content-between h-100">
            <p class="description flex-grow-1 mb-3">{component.meta.description}</p>
            {#if component.meta.keywords}
                <p class="topics mb-0 ms-n2">
                    {#each component.meta.keywords as keyword}
                        <span class={`badge fw-normal bg-body-tertiary text-success me-2`}>{keyword}</span>
                    {/each}
                </p>
            {/if}
            {#if component.pipelines}
                <div class="text-body-secondary align-self-bottom included-in">
                    <span class="text-small">Included in:</span>
                    {#if collapsePipelines}
                        {#each component.pipelines.slice(0, 3) as pipeline}
                            <a
                                class="badge fw-normal border border-warning-subtle bg-warning-subtle text-body me-2 text-decoration-none"
                                href={"/" + pipeline.name + "/" + pipeline.version + "/"}>{pipeline.name}</a
                            >
                        {/each}
                        {#if component.pipelines.length > 3}
                            <span
                                class="text-small cursor-pointer"
                                title="click to show all pipelines"
                                onclick={() => (collapsePipelines = !collapsePipelines)}
                                onkeydown={() => (collapsePipelines = !collapsePipelines)}
                                tabindex="0"
                                role="button"
                                data-bs-toggle="tooltip"
                                data-bs-delay="500">+{component.pipelines.length - 3} more pipelines</span
                            >
                        {/if}
                    {:else}
                        {#each component.pipelines as pipeline}
                            <a
                                class="badge fw-normal border border-warning-subtle bg-warning-subtle text-body me-2"
                                href={"/" + pipeline.name + "/" + pipeline.version + "/"}>{pipeline.name}</a
                            >
                        {/each}
                    {/if}
                </div>
            {/if}
            {#if component.type !== "module" && component.meta.components}
                <div class="text-body-secondary align-self-bottom components">
                    <span class="text-small">Includes:</span>
                    {#each component.meta.components as sub_component}
                        <span
                            class={`badge fw-normal border border-info-subtle bg-info-subtle text-body me-2 ${component.type}-topic`}
                            >{sub_component}</span
                        >
                    {/each}
                </div>
            {/if}
            {#if component.subworkflows}
                <span class="text-body-secondary align-self-bottom">
                    <span class="text-small">Part of:</span>
                    {#each component.subworkflows as subworkflow}
                        <a
                            class="badge fw-normal border border-info-subtle bg-info-subtle text-body me-2"
                            href={"/subworkflows/" + subworkflow + "/"}>{subworkflow}</a
                        >
                    {/each}
                </span>
            {/if}
            {#if component.meta.modules}
                <span class="text-body-secondary align-self-bottom"
                    >Includes:
                    {#each component.meta.modules as module}
                        <a
                            class="badge fw-normal border border-info-subtle bg-info-subtle text-body me-2"
                            href={"/modules/" + module.replace("/", "_") + "/"}>{module}</a
                        >
                    {/each}
                </span>
            {/if}
        </div>
    {/snippet}
</ListingCard>

<style>
    .card {
        max-width: 30rem;
    }
    .card-header {
        background-color: transparent;
        word-break: break-word;
    }
    .included-in .badge:last-child {
        margin-right: 0;
    }
</style>
