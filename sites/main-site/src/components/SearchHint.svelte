<script lang="ts">
    import { SearchQuery } from "@components/store";
    import type { SearchExample } from "@utils/search";

    let { examples = [] }: { examples: SearchExample[] } = $props();
</script>

<div class="text-center">
    <div class="text-secondary-emphasis px-3 mb-3 text-small">
        <strong>Hint:</strong> Scope your search to one field with a prefix:
        {#each examples as example, i (example.query)}
            <button
                type="button"
                class="fw-normal border-0 bg-body-tertiary text-secondary-emphasis font-monospace"
                onclick={() => SearchQuery.set(example.query)}>{example.query}</button
            >({example.label}){i < examples.length - 1 ? "," : "."}
        {/each}
    </div>
</div>
