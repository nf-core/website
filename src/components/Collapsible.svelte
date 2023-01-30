<script>
  import {showHelp} from './store.js'
  export let buttonText = "help text";
  export let buttonClass = "btn-outline-info me-auto";
  export let textClass = "";
  $: show = $showHelp;
</script>

<div>
  <div class="">
    <button class={"btn " + buttonClass} class:open={show} type="button" on:click={() => (show = !show)}>
      {buttonText}
    </button>
  </div>
  <div class={"collapse " + textClass} class:show>
    <slot />
  </div>
</div>

<style lang="scss">
  .btn.open {
    border-bottom-left-radius: 0;
    border-bottom-right-radius: 0;
  }
</style>
