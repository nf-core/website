<script>
    import SchemaListingGroup from '@components/schema/SchemaListingGroup.svelte';
    import { onMount } from 'svelte';
    import { currentHeading } from '@components/store';

    export let schema;
    onMount(() => {
        const observer = new IntersectionObserver(
            (entries) => {
                entries.forEach((entry) => {
                    if (entry.isIntersecting) {
                        currentHeading.set(entry.target.id);
                    }
                });
            },
            {
                rootMargin: '0px 0px -92% 0px',
            }
        );
        Object.entries(schema.definitions).forEach((heading) => {
            const element = document.querySelector('#' + heading[0].replaceAll('_', '-'));
            console.log(heading);
            observer.observe(element);
        });
    });
</script>

<div class="schema-listing">
    {#if Object.entries(schema.definitions).length > 0}
        <div class="d-flex flex-column">
            {#each Object.entries(schema.definitions) as [id, definition] (id)}
                <SchemaListingGroup {definition} {id} />
            {/each}
        </div>
    {:else}
        <div class="alert alert-warning mt-3" role="alert">
            <h4 class="text-warning">No nextflow_schema.json file found!</h4>
            <p>
                It seems like there is no nextflow_schema.json file with parameters defined for this version of the
                pipeline. Try a newer version.
            </p>
        </div>
    {/if}
</div>

<style>
    .alert-warning {
        max-width: 50rem;
    }
</style>
