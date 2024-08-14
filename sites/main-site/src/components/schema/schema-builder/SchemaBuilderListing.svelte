<script lang="ts">
    import { onMount } from 'svelte';

    // import { dnd } from 'svelte-pragmatic-list';
    import SchemaBuilderListingElement from './SchemaBuilderListingElement.svelte';
    import SchemaBuilderToolbar from './SchemaBuilderToolbar.svelte';
    let schema = {};

    let properties = {};
    let id = '';

    // rename key in the schema
    const updateSchemaKey = (e, oldKey) => {
        // handle edge cases
        const newKey = e.target.value;
        if (newKey === oldKey) return;
        if (newKey in properties) return;
        if (newKey.length == 0) return;
        const newSchema = {};
        const keys = Object.keys(properties);

        // rename the key in the schema, but keep the order the sam
        for (let i = 0; i < keys.length; i++) {
            if (keys[i] === oldKey) {
                newSchema[newKey] = properties[oldKey];
            } else {
                newSchema[keys[i]] = properties[keys[i]];
            }
        }
        schema.definitions['input_output_options'].properties = newSchema;
    };
    // post schema to the server
    const postSchema = async () => {
        const response = await fetch(`http://localhost:8000/process-schema?id=${id}`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                message: 'Schema created successfully',
            },
            body: JSON.stringify({ status: 'web_builder_edited', schema }),
        });
        console.log(await response.json());
    };
    let loading = true;
    onMount(async () => {
        // read id from URL
        const url = new URL(window.location.href);
        id = url.searchParams.get('id') || '';
        if (!id) {
            console.error('No id found in the URL');
            return;
        }
        try {
            const response = await fetch(`http://localhost:8000/process-schema?id=${id}`);
            const data = await response.json();

            schema = data.data.schema;
            console.log('schema', schema);

            properties = schema.definitions['input_output_options'].properties;
            loading = false;
        } catch (error) {
            console.error('Error fetching schema:', error);
            loading = false;
            // Handle the error, maybe show an error message to the user
        }
    });
    console.log('properties', schema);
</script>

<SchemaBuilderToolbar on:update={({ detail }) => (schema = detail)}>
    <button id="finish" class="btn btn-primary" on:click={postSchema}>Finish</button>
</SchemaBuilderToolbar>

<div class="border rounded-md p-2">
    {#if loading}
        <p>Loading schema...</p>
    {:else}
        {#each Object.keys(properties) as key (key)}
            <SchemaBuilderListingElement bind:schema={properties} name={key}>
                <label
                    ><span class:d-none={key.length == 0}>ID</span>
                    <input
                        type="text"
                        class="font-monospace"
                        value={key}
                        placeholder="ID"
                        on:change={(e) => {
                            updateSchemaKey(e, key);
                        }}
                    />
                </label>
                {key}
            </SchemaBuilderListingElement>
        {/each}
        <div data-dnd-indicator hidden style="height: 20px; background-color: red; width: 100%;"></div>
    {/if}
</div>
<pre class="border rounded p-2">
  {JSON.stringify(schema, null, 2)}
</pre>
