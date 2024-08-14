<script lang="ts">
    export let schema = {};

    const addField = () => {
        console.log('add field');
        schema.definitions['new-group'] = {
            properties: {
                'new-field': {
                    type: 'string',
                    description: '',
                    fa_icon: '',
                    help_text: '',
                    pattern: '',
                },
            },
        };
    };
</script>

<div class="d-flex justify-content-between border rounded p-2 mb-2">
    <div class="d-flex">
        <button class="btn btn-outline-secondary me-2" on:click={addField}>Add Field</button>
        <button class="btn btn-outline-secondary me-2 disabled">Add Group</button>
    </div>

    <slot />
</div>
