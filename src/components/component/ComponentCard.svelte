<script lang="ts">
    import ListingCard from '@components/ListingCard.svelte';
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
    const href = '/' + component.type + 's/' + component.name + '/';
</script>

<ListingCard>
    <div slot="card-header">
        <a {href}>{@html component.name.replace('_', '_<wbr>')} </a>
        <small class="gh-stats text-small" />
    </div>
    <div slot="card-body" class="d-flex flex-column justify-content-between h-100">
        {#if component.meta.keywords}
            <p class="topics">
                {#each component.meta.keywords as keyword}
                    <span class="badge fw-normal bg-body-tertiary text-success me-2">{keyword}</span>
                {/each}
            </p>
        {/if}

        <p class="description flex-grow-1">{component.meta.description}</p>
        {#if component.pipelines}
            <span class="text-body-secondary align-self-bottom"
                >Included in:
                {#each component.pipelines as pipeline}
                    <a
                        class="badge fw-normal border border-warn-subtleing bg-warning-subtle text-body me-2"
                        href={'/' + pipeline.name + '/' + pipeline.version + '/'}>{pipeline.name}</a
                    >
                {/each}
            </span>
        {/if}
        {#if component.subworkflows}
            <span class="text-body-secondary align-self-bottom"
                >Part of:
                {#each component.subworkflows as subworkflow}
                    <a
                        class="badge fw-normal border border-info-subtle bg-info-subtle text-body me-2"
                        href={'/subworkflows/' + subworkflow + '/'}>{subworkflow}</a
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
                        href={'/modules/' + module.replace('/', '_') + '/'}>{module}</a
                    >
                {/each}
            </span>
        {/if}
    </div>
</ListingCard>

<style>
    .card {
        max-width: 30rem;
    }
    .card-header {
        background-color: transparent;
        word-break: break-word;
    }
</style>
