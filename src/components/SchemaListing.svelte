<script>
  import SchemaListingGroup from "./SchemaListingGroup.svelte";

  export let schema;
</script>

<div class="schema-listing">
  {#if schema.definitions}
    <div class="d-flex flex-column">
      {#each Object.entries(schema.definitions) as [id, definition] (id)}
        <SchemaListingGroup {definition} {id}/>
      {/each}
    </div>
  {/if}
</div>
