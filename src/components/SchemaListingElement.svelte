<script>
  import {showHidden} from './store.js'

  export let title;
  export let property;
  const id = title.replace(" ","-")

  import { marked } from "marked";
  import Collapsible from "./Collapsible.svelte";
  let help_text = "";
  if (property.help_text) {
    help_text = marked.parse(property.help_text);
  }
</script>

<div class="property row border border-secondary py-2 mx-md-2 justify-content-between" class:collapse={property.hidden} class:show={$showHidden}>
  <div class="col-12 col-md-3 title border-right border-secondary overflow-x-scroll text-nowrap">
    {#if property.fa_icon}
      <i class="fa fa-fw me-2 {property.fa_icon}" />
    {/if}
    <code>--{title}</code>
  </div>
  <div class="col description ">
    {property.description}
  </div>
  <div class="col-12 col-md-3 text-nowrap d-flex flex-column align-items-end justify-content-between">
    {#if property.hidden}
      <span class="badge text-bg-warning">hidden</span>
    {/if}
    <div class="text-muted">
      type: <code>{property.type}</code>
    </div>
    {#if property.enum}
      <select value={property.default}>
        {#each property.enum as value}
          <option {value}>{value + (value === property.default ? " (default)" : "")}</option>
        {/each}
      </select>
    {/if}
    {#if property.default && !property.enum}
      <div class="default w-100 text-end text-muted overflow-x-scroll">
        default: <code>{property.default}</code>
      </div>
    {/if}
  </div>
  {#if property.help_text}
    <div class="row my-2 mx-1 w-100 ms-auto px-0">
      <Collapsible>

      <div {id} class="p-2 text-info-emphasis bg-info-subtle border border-info-subtle rounded-3">
        {@html help_text}
      </div>
      </Collapsible>
    </div>
  {/if}
</div>

<style lang="scss">

  .rounded-3{
    border-top-left-radius: 0 !important;
  }
</style>
