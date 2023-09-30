<script>
    import SchemaListingElement from '@components/schema/SchemaListingElement.svelte';
    import { showHidden } from '@components/store';

    export let definition;
    export let id;

    //go through definition.properities and check if all hidden are set to true
    let hidden = true;
    for (const [title, property] of Object.entries(definition.properties)) {
        if (!property.hidden) {
            hidden = false;
            break;
        }
    }
    //get definition.required and assign required=true to the corresponding defintion.properties

    let required = definition.required;
    if (!required) required = [];
    required.forEach((item) => {
        if (definition.properties[item]) {
            definition.properties[item].required = true;
        }
    });
</script>

<div class="card my-2" class:collapse={hidden} class:show={$showHidden}>
    <div class="card-header position-sticky bg-body-secondary">
        <h2 class="card-title text-success scroll-target my-1" id={id.replaceAll('_', '-')}>
            <a
                class="text-decoration-none text-success"
                aria-hidden="true"
                tabindex="-1"
                href={'#' + id.replaceAll('_', '-')}
                ><i class="ms-1 fas invisible" aria-hidden="true" />{#if definition.fa_icon}
                    <i class="fa fa-fw me-2 {definition.fa_icon}" />
                {/if}
                {definition.title}
            </a>
        </h2>
    </div>
    <div class="card-body pb-0">
        <p class="mb-0">{definition.description ? definition.description : ''}</p>
        {#if definition.properties}
            <div class="properties">
                {#each Object.entries(definition.properties) as [title, property] (title)}
                    <SchemaListingElement {title} {property} />
                {/each}
            </div>
        {/if}
    </div>
</div>

<style lang="scss">
    .properties > :global(.row) {
        margin-top: -1px; // avoid doubled borders
    }
</style>
