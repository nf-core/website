<script lang="ts">
    import sanitizeHtml from 'sanitize-html';

    export let name: string = '';

    export let schema: any = {};
    export let required: boolean = false;

    $: schemaEntry = schema[name];

    const updateSchema = (e, propName) => {
        schema[name][propName] = e.target.value;
    };
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
            <input
                type="text"
                class=""
                data-param_key="description"
                value={sanitizeHtml(schemaEntry.description)}
                placeholder="description"
                on:input={(e) => updateSchema(e, 'description')}
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
            Required
        </label>
    </div>
</div>
