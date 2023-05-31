<script lang="ts">
    export let component: {
        name: string;
        path: string;
        type: string;
        meta: {
            description: string;
            name: string;
            keywords?: string[];
            modules?: string[];
        };
        pipelines?: {
            name: string;
            version: string;
        }[];
        subworkflows?: string[];
    };
    const href = '/' + component.type + 's/' + component.name;
</script>

<div class="card flex-fill m-2">
    <div class="card-header border-bottom-0 pb-0">
        <h2 class="mb-0 d-flex justify-content-between align-items-center">
            <a {href}>{@html component.name.replace('_', '_<wbr>')} </a>
            <small class="gh-stats text-small" />
        </h2>
    </div>
    <div class="card-body pt-0 d-flex flex-column">
        {#if component.meta.keywords}
            <p class="topics">
                {#each component.meta.keywords as keyword}
                    <span class="badge bg-body-tertiary text-success me-2">{keyword}</span>
                {/each}
            </p>
        {/if}

        <p class="description flex-grow-1">{component.meta.description}</p>
        {#if component.pipelines}
            <span class="text-body-secondary align-self-bottom"
                >Included in:
                {#each component.pipelines as pipeline}
                    <a class="badge text-bg-warning me-2" href={'/' + pipeline.name + '/' + pipeline.version}
                        >{pipeline.name}</a
                    >
                {/each}
            </span>
        {/if}
        {#if component.subworkflows}
            <span class="text-body-secondary align-self-bottom"
                >Part of:
                {#each component.subworkflows as subworkflow}
                    <a class="badge text-bg-info me-2" href={'/subworkflows/' + subworkflow}>{subworkflow}</a>
                {/each}
            </span>
        {/if}
        {#if component.meta.modules}
            <span class="text-body-secondary align-self-bottom"
                >Contains:
                {#each component.meta.modules as module}
                    <a class="badge text-bg-info me-2" href={'/modules/' + module.replace('/', '_')}>{module}</a>
                {/each}
            </span>
        {/if}
    </div>
</div>

<style>
    .card {
        max-width: 30rem;
    }
    .card-header {
        background-color: transparent;
        word-break: break-word;
    }
</style>
