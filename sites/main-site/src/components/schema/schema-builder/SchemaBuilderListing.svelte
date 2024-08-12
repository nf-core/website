<script lang="ts">
    import { dnd } from 'svelte-pragmatic-list';
    import SchemaBuilderListingElement from './SchemaBuilderListingElement.svelte';
    import SchemaBuilderToolbar from './SchemaBuilderToolbar.svelte';

    export let schema = {};
    export let id: string = '';
    let items =
        schema && schema.schema && schema.schema['definitions'] && Object.values(schema.schema['definitions'])
            ? Object.entries(Object.values(schema.schema['definitions'])[1]['properties'])
            : [];
</script>

<SchemaBuilderToolbar {id} {schema} />
<div
    use:dnd={{
        items,
        onChange: (newItems) => (items = newItems),
    }}
    class="border rounded-md p-2"
>
    {#each items as item}
        <SchemaBuilderListingElement bind:schemaEntry={item} />
        {item[0]}
    {/each}
    <div data-dnd-indicator hidden style="height: 20px; background-color: red; width: 100%; "></div>
</div>
