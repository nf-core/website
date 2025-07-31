<script>
    import SchemaListingElement from "@components/schema/SchemaListingElement.svelte";
    import { showHidden } from "@components/store";
    import Markdown from "@components/markdown/Markdown.svelte";

    let { definition = $bindable(), id } = $props();

    let processedProperties = $derived(() => {
        if (!definition.properties) return {};

        const properties = { ...definition.properties };
        const required = definition.required || [];

        // Mark required properties without mutating the original
        required.forEach((item) => {
            if (properties[item]) {
                properties[item] = { ...properties[item], required: true };
            }
        });

        return properties;
    });

    // Calculate if all properties are hidden using $derived
    let hidden = $derived(
        !definition.properties || Object.values(definition.properties).every((property) => property.hidden),
    );
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
                {@html definition.description_rendered}
            {/if}
        </p>
        {#if processedProperties()}
            <div class="properties">
                {#each Object.entries(processedProperties()) as [title, property] (title)}
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
