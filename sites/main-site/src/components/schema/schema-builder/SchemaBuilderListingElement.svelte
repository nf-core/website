<script lang="ts">
    import sanitizeHtml from 'sanitize-html';

    export let name: string = '';

    export let schema: any = {};
    export let required: boolean = false;

    $: schemaEntry = schema[name];

    const updateSchema = (e, propName) => {
        schema[name][propName] = e.target.value;
    };

    const schemaTypes = ['string', 'number', 'integer', 'boolean'];
</script>

<div class="row border">
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
        <label
            >Description
            <textarea
                class=""
                data-param_key="description"
                value={sanitizeHtml(schemaEntry.description)}
                placeholder="description"
                on:input={(e) => updateSchema(e, 'description')}
            />
        </label>
    </div>
    <div class="col">
        <label
            >Help Text
            <textarea
                class=""
                data-param_key="help_text"
                value={sanitizeHtml(schemaEntry.help_text)}
                placeholder="help text"
                on:input={(e) => updateSchema(e, 'help_text')}
            />
        </label>
    </div>
    <div class="col-auto">
        <label
            >Type
            <select class="param_key param_type" data-param_key="type" on:change={(e) => updateSchema(e, 'type')}>
                {#each schemaTypes as type}
                    <option selected={type == schemaEntry.type} value={type}>{type}</option>
                {/each}
            </select>
        </label>
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
