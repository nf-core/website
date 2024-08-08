<script lang="ts">
    import { dnd } from 'svelte-pragmatic-list';
    import SchemaBuilderListingElement from './SchemaBuilderListingElement.svelte';
    export let schema: {
        $schema: string;
        $id: string;
        title: string;
        description: string;
        type: string;
        definitions: {
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
        definitions: {},
    };
    console.log('schema', schema);
    let items =
        schema && schema['definitions'] && Object.values(schema['definitions'])
            ? Object.entries(Object.values(schema['definitions'])[0]['properties'])
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
