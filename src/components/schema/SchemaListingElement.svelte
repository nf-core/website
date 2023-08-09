<script>
    import { showHidden } from '@components/store';
    import Collapsible from '@components/Collapsible.svelte';
    import Markdown from '@components/markdown/Markdown.svelte';

    export let title;
    export let property;
    const id = title.replace(' ', '-');
</script>

<div
    class="property row border-bottom py-3 mx-md-1 justify-content-between"
    class:collapse={property.hidden}
    class:show={$showHidden}
>
    <div id={title} class="col-12 col-md-3 title border-right border-secondary text-nowrap p-0 overflow-x-scroll">
        <a class="text-decoration-none" aria-hidden="true" tabindex="-1" href={'#' + title}
            ><i class="ms-1 fas invisible" aria-hidden="true" />
            <span class="">
                {#if property.fa_icon}
                    <i class="fa fa-fw {property.fa_icon}" />
                {/if}
                <code>--{title}</code>
            </span>
        </a>
    </div>
    <div class="col description text-small">
        <Markdown md={property.description} />
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
                    <option {value}>{value + (value === property.default ? ' (default)' : '')}</option>
                {/each}
            </select>
        {/if}
        {#if property.default && !property.enum}
            <div class="default w-100 text-end text-body-secondary overflow-x-scroll">
                default: <code>{property.default}</code>
            </div>
        {/if}
    </div>
    {#if property.help_text}
        <div class="row d-flex mt-2 mx-0 w-100 px-0 gx-3 gx-md-4">
            <Collapsible>
                <div {id} class="p-2 px-3 text-body bg-secondary-subtle border border-secondary rounded-3">
                    <Markdown md={property.help_text} />
                </div>
            </Collapsible>
        </div>
    {/if}
</div>

<style lang="scss">
    .rounded-3 {
        border-top-right-radius: 0 !important;
        margin-top: -1pt; // avoid doubled borders
    }
    .property {
        margin-bottom: -1pt;
    }
</style>
