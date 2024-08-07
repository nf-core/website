<script lang="ts">
    import { dnd } from 'svelte-pragmatic-list';
    import SchemaBuilderListingElement from './SchemaBuilderListingElement.svelte';
    export let schema: {
        $schema: string;
        $id: string;
        title: string;
        description: string;
        type: string;
        $defs: {
            [key: string]: {
                type: string;
                properties: {
                    [key: string]: {
                        id: string;
                        description: string;
                        type: string;
                        default: string;
                        required: boolean;
                        hidden: boolean;
                        fa_icon: string;
                        help_text: string;
                    };
                };
            };
        };
    } = {
        $schema: '',
        $id: '',
        title: '',
        description: '',
        type: '',
        $defs: {},
    };

    let items =
        schema && schema['$defs'] && Object.values(schema['$defs'])
            ? Object.entries(Object.values(schema['$defs'])[0]['properties'])
            : [];
    console.log('items', items);
</script>

<div
    use:dnd={{
        items,
        onChange: (newItems) => (items = newItems),
    }}
    class="border rounded-md p-2"
>
    {#each items as item}
        <SchemaBuilderListingElement schemaEntry={item} />
    {/each}
    <div data-dnd-indicator hidden style="height: 20px; background-color: red; width: 100%; "></div>
</div>
