<script lang="ts">
    import { showHelp, showHidden } from '@components/store';
    import SchemaBuilderListingElement from './SchemaBuilderListingElement.svelte';
    import SchemaListingElement from '@components/schema/SchemaListingElement.svelte';
    export let name: string = '';
    export let schema: any = {};

    showHelp.set(true);
    showHidden.set(true);

    $: group = schema.definitions[name];

    let fieldAdded = false;
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
    const addField = () => {
        let newName = 'new_field';
        // check if field already exists and add number if it does (or add 1 to the number if also a number already exists)
        let i = 1;
        while (newName in schema.definitions[name].properties) {
            if (i === 1) {
                newName = `${newName}_${i}`;
            } else {
                newName = newName.replace(/_\d+$/, `_${i}`);
            }
            i++;
        }

        schema.definitions[name].properties = {
            [newName]: {
                type: 'string',
                description: 'New Field',
            },
            ...schema.definitions[name].properties,
        };
        fieldAdded = true;
        // reset fieldAdded after 1 second
        setTimeout(() => {
            fieldAdded = false;
        }, 1000);
    };

    const removeField = (property) => {
        let result = confirm('Are you sure you want to remove this field?');
        if (!result) return;
        delete schema.definitions[name].properties[property];
        schema.definitions[name].properties = { ...schema.definitions[name].properties };
    };
</script>

<div class="border border-success bg-bo">
    <div class="d-flex justify-content-between border rounded p-2 mb-2">
        <h3><i class={`fa ${group.fa_icon} me-2 fa-fw`}></i>{name}</h3>
        <slot name="delete-group"></slot>
        <button class="btn btn-outline-secondary" on:click={addField}>
            {#if fieldAdded}
                <span><i class="fas fa-check"></i> Field Added</span>
            {:else}
                <span><i class="fas fa-plus"></i> Add Field</span>
            {/if}
        </button>
    </div>
    {#if group.description}
        <p>{group.description}</p>
    {/if}
    {#each Object.keys(group.properties) as key (key)}
        <div class="grid justify-items-stretch">
            <div class="g-col-12 g-col-md-6">
                <SchemaListingElement title={key} property={group.properties[key]} />
            </div>
            <div class="g-col-12 g-col-md-6 bg-body-tertiary">
                <SchemaBuilderListingElement
                    bind:schema={group.properties}
                    name={key}
                    required={group.required?.includes(key)}
                >
                    <div class="form-floating">
                        <input
                            id="param_id"
                            type="text"
                            class="font-monospace form-control"
                            value={key}
                            placeholder="ID"
                            on:change={(e) => {
                                updateSchemaKey(e, key);
                            }}
                        />
                        <label for="param_id">ID </label>
                    </div>
                    <div class="col-auto ms-auto" slot="delete">
                        <button class="btn btn-outline-danger" on:click={() => removeField(key)}>
                            <i class="fas fa-trash"></i>
                        </button>
                    </div>
                </SchemaBuilderListingElement>
            </div>
        </div>
    {/each}
</div>
