<script>
    import SchemaListingElement from "@components/schema/SchemaListingElement.svelte";
    import { showHidden } from "@components/store";
    import Markdown from "@components/markdown/Markdown.svelte";

    let { definition = $bindable(), id } = $props();

    // Go through definition.properties and check if all hidden are set to true
    let hidden = $state(true);

    $effect(() => {
        // Reset hidden state when definition changes
        hidden = true;
        console.log(definition);
        if (definition.properties) {
            for (const [title, property] of Object.entries(definition.properties)) {
                if (!property.hidden) {
                    hidden = false;
                    break;
                }
            }
        }

        // Get definition.required and assign required=true to the corresponding definition.properties
        let required = definition.required || [];
        required.forEach((item) => {
            if (definition.properties && definition.properties[item]) {
                definition.properties[item].required = true;
            }
        });
    });
</script>

<div class="card my-2" class:collapse={hidden} class:show={$showHidden}>
    <div class="card-header position-sticky bg-body-secondary">
        <h2 class="card-title text-success scroll-target my-1" id={id?.replaceAll("_", "-")}>
            <a
                class="text-decoration-none text-success"
                aria-hidden="true"
                tabindex="-1"
                href={"#" + id?.replaceAll("_", "-")}
                ><i class="ms-1 fas fa-xs invisible" aria-hidden="true"></i>{#if definition.fa_icon}
                    <i class="fa fa-fw me-2 {definition.fa_icon}"></i>
                {/if}
                {definition.title}
            </a>
        </h2>
    </div>
    <div class="card-body pb-0">
        <p class="mb-0">
            {#if definition.description}
                <Markdown md={definition.description} />
            {/if}
        </p>
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
