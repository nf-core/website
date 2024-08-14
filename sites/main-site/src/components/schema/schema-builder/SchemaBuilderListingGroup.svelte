<script lang="ts">
    import SchemaBuilderListingElement from './SchemaBuilderListingElement.svelte';

    export let name: string = '';
    export let schema: any = {};

    $: group = schema.definitions[name];
    const updateSchemaKey = (e, oldKey) => {
        // handle edge cases
        const newKey = e.target.value;
        if (newKey === oldKey) return;
        if (newKey in schema.definitions[name].properties) return;
        if (newKey.length == 0) return;
        const newSchema = {};
        const keys = Object.keys(schema.definitions[name].properties);

        // rename the key in the schema, but keep the order the sam
        for (let i = 0; i < keys.length; i++) {
            if (keys[i] === oldKey) {
                newSchema[newKey] = schema.definitions[name].properties[oldKey];
            } else {
                newSchema[keys[i]] = schema.definitions[name].properties[keys[i]];
            }
        }
        schema.definitions[name].properties = newSchema;
    };
</script>

<div class="border border-success">
    <h3><i class={`fa ${group.fa_icon} me-2 fa-fw`}></i>{name}</h3>
    {#if group.description}
        <p>{group.description}</p>
    {/if}
    {#each Object.keys(group.properties) as key (key)}
        <SchemaBuilderListingElement bind:schema={group.properties} name={key} required={group.required?.includes(key)}>
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
</div>
