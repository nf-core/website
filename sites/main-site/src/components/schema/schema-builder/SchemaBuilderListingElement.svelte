<script lang="ts">
    import sanitizeHtml from 'sanitize-html';

    export let name: string = '';

    export let schema: any = {};
    export let required: boolean = false;

    $: schemaEntry = schema[name];

    const updateSchema = (event, propName: string) => {
        schema[name][propName] = event.target.value;
    };

    const schemaTypes = ['string', 'number', 'integer', 'boolean'];
</script>

<div class="row border align-items-center pt-2">
    <div data-dnd-handle class="col-auto align-self-center border-end">
        <i class="fas fa-grip-vertical"></i>
    </div>
    <button class="btn col-auto align-self-center">
        <i class={schemaEntry.fa_icon}></i>
    </button>
    <div class="col">
        <slot />
    </div>
    <div class="d-sm-none w-100"></div>
    <div class="col">
        <div class="form-floating">
            <textarea
                class="description-textarea form-control"
                id="description"
                placeholder="description"
                data-param_key="description"
                on:change={(e) => updateSchema(e, 'description')}
                on:keydown={(e) => updateSchema(e, 'description')}>{sanitizeHtml(schemaEntry.description)}</textarea
            >
            <label for="description">Description</label>
        </div>
    </div>

    <div class="col">
        <div class="form-floating">
            <textarea
                class="help_text-textarea form-control"
                id="help_text"
                placeholder="help_text"
                data-param_key="help_text"
                on:change={(e) => updateSchema(e, 'help_text')}
                on:keydown={(e) => updateSchema(e, 'help_text')}>{sanitizeHtml(schemaEntry.help_text)}</textarea
            >
            <label for="help_text">help text</label>
        </div>
    </div>
    <div class="col-auto">
        <label for="type">Type</label>

        <select id="type" class="param_key param_type" data-param_key="type" on:change={(e) => updateSchema(e, 'type')}>
            {#each schemaTypes as type}
                <option selected={type == schemaEntry.type} value={type}>{type}</option>
            {/each}
        </select>
    </div>
    <div class="d-sm-none w-100"></div>

    <div class="col">
        <label
            >Default
            <input
                type="text"
                class=""
                data-param_key="default"
                value={schemaEntry.default ?? ''}
                placeholder="default"
                on:input={(e) => updateSchema(e, 'default')}
            />
        </label>
    </div>
    <!-- check box for required -->
    <div class="col-auto align-self-center">
        <label>
            <input
                type="checkbox"
                class="form-check-input"
                data-param_key="required"
                bind:checked={required}
                on:change={(e) => updateSchema(e, 'required')}
            />
            required
        </label>
    </div>
    <div class="col-auto align-self-center">
        <label>
            <input
                type="checkbox"
                class="form-check-input"
                data-param_key="required"
                bind:checked={schemaEntry.hidden}
                on:change={(e) => updateSchema(e, 'hidden')}
            />
            hidden
        </label>
    </div>
</div>

<style lang="scss">
    .description-textarea,
    .help_text-textarea {
        height: 5rem;
    }
</style>
