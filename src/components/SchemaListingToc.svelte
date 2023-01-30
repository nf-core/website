
<script>
  import {showHidden,showHelp} from './store.js'

  export let schema;
  // console.log(schema.definitions);
</script>

<div class="toc">
  <h2>On this page</h2>
  <button type="button" class="btn btn-outline-secondary" on:click={()=>showHidden.set(!$showHidden)}>Show hidden</button>
  <button type="button" class="btn btn-outline-secondary" on:click={()=>showHelp.set(!$showHelp)}>Show all help</button>
  
  <!-- <p  class="text-success">{schema.description}</p>
  {#if schema.definitions}
    <div class="d-flex flex-column">

      {#each Object.entries(schema.definitions) as [_, definition]}
        
      {/each}
    </div>
  {/if} -->
</div>
