<script>
    export let module;

    const href = '/modules/' + module.name;
</script>

<div class="card flex-fill m-2">
    <div class="card-header border-bottom-0 pb-0">
        <h2 class="mb-0 d-flex justify-content-between align-items-center">
            <a {href}>{@html module.name.replace('_', '_<wbr>')} </a>
            <small class="gh-stats text-small" />
        </h2>
    </div>
    <div class="card-body pt-0 d-flex flex-column">
        {#if module.meta.keywords}
            <p class="topics">
                {#each module.meta.keywords as keyword}
                    <span class="badge bg-body-tertiary text-success me-2">{keyword}</span>
                {/each}
            </p>
        {/if}

        <p class="description flex-grow-1">{module.meta.description}</p>
        {#if module.pipelines}
            <span class="text-muted align align-self-bottom"
                >Used in:
                {#each module.pipelines as pipeline}
                    <a class="badge text-bg-warning me-2" href={'/' + pipeline.name + '/' + pipeline.version}
                        >{pipeline.name}</a
                    >
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
