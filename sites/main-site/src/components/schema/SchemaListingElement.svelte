<script lang="ts">
    import { showHidden } from "@components/store";
    import Collapsible from "@components/Collapsible.svelte";

    let { title, property } = $props();
    const id = $derived(title.replace(" ", "-"));

    let longPattern = $derived(
        (() => {
            const pattern = property?.pattern;

            // explicitly handle patterns which are an enum work around, i.e. they have multiple values, eg. "^(foo|bar)$"
            if (
                pattern &&
                pattern.startsWith("^(") &&
                (pattern.endsWith(")$") || pattern.endsWith(")*$")) &&
                pattern.includes("|")
            ) {
                const matches = pattern.match(/\b(\w+)\b/g);
                return matches || [];
            }
            return [];
        })(),
    );

    let formattedLongPattern = $derived(
        longPattern.length ? "<code>" + longPattern.join("</code>, <code>") + "</code>" : "",
    );
</script>

<div
    class="property row border-bottom py-3 mx-md-1 justify-content-between align-items-center"
    class:collapse={property.hidden}
    class:show={$showHidden}
>
    <div id={title} class="col-12 col-md-3 title border-right border-secondary text-nowrap p-0 pe-2">
        <a class="text-decoration-none d-block overflow-x-scroll" aria-hidden="true" tabindex="-1" href={"#" + title}
            ><i class="ms-1 fas invisible" aria-hidden="true"></i>
            <span class="">
                {#if property.fa_icon}
                    <i class="fa fa-fw {property.fa_icon}"></i>
                {/if}
                <code>--{title}</code>
            </span>
        </a>
    </div>
    <div class="col description text-small">
        {@html property.description_rendered}
    </div>
    <div class="col-12 col-md-3 text-nowrap d-flex flex-column align-items-end justify-content-between">
        {#if property.hidden}
            <span class="badge warning border border-warning-subtle fw-normal bg-warning-subtle mb-1">hidden</span>
        {/if}
        {#if property.required}
            <span class="badge text-bg-warning mb-1">required</span>
        {/if}
        <div class="text-body-secondary">
            type: <code>{property.type}</code>
        </div>
        {#if property.enum}
            <select class="form-select mt-2" value={property.default} aria-label="Parameter enum options">
                {#each property.enum as value}
                    <option {value}>{value + (value === property.default ? " (default)" : "")}</option>
                {/each}
            </select>
        {/if}
        {#if property.default && !property.enum}
            <div class="default w-100 text-end text-body-secondary overflow-x-scroll">
                default: <code>{property.default}</code>
            </div>
        {/if}
        {#if property.pattern && !longPattern.length}
            <div class="default w-100 text-end text-body-secondary overflow-x-scroll">
                pattern: <code>{property.pattern}</code>
            </div>
        {/if}
    </div>
    {#if property.help_text}
        <div class="row d-flex mt-2 mx-0 w-100 px-0 gx-3 gx-md-4 help-text">
            <Collapsible>
                <div {id} class="p-2 px-3 border border-secondary rounded-3">
                    {@html property.help_text_rendered}
                    {#if longPattern.length}
                        This parameter must be a combination of the following values:
                        {@html formattedLongPattern}
                    {/if}
                </div>
            </Collapsible>
        </div>
    {/if}
</div>

<style lang="scss">
    .help-text .rounded-3 {
        border-top-left-radius: 0 !important;
        margin-top: -0.75pt; // avoid doubled borders
    }
    .property {
        margin-bottom: -1pt;
    }

    .help-text :global(pre code) {
        padding-left: 1rem;
    }
</style>
