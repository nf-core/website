<script lang="ts">
    import sanitizeHtml from 'sanitize-html';

    export let name: string = '';

    export let schema: any = {};
    export let required: boolean = false;

    $: schemaEntry = schema[name];

    const updateSchema = (event, propName: string) => {
        let value = event.target.value;
        if (propName === 'enum') {
            value = value.split(',').map((item) => item.trim());
        }
        schema[name][propName] = value;
        if (schemaEntry.pattern && (propName === 'default' || propName === 'enum')) {
            checkPattern(propName);
        }
        if (
            schemaEntry.multipleOf &&
            (propName === 'default' || propName === 'minimum' || propName === 'maximum' || propName === 'multipleOf')
        ) {
            checkMultipleOf();
        }
    };
    const schemaTypes = ['string', 'number', 'integer', 'boolean'];
    let validation: {} = {};
    const checkPattern = (key) => {
        // check if default value and enum values match the pattern
        if (schemaEntry.pattern) {
            if (key === 'default' && schemaEntry.default) {
                if (!schemaEntry.default.match(schemaEntry.pattern)) {
                    validation['default'] = false;
                } else {
                    validation['default'] = true;
                }
            } else if (key === 'enum' && schemaEntry.enum) {
                for (let i = 0; i < schemaEntry.enum.length; i++) {
                    if (!schemaEntry.enum[i].match(schemaEntry.pattern)) {
                        validation['enum'] = false;
                        break;
                    }
                }
                validation['enum'] = true;
            }
        }
    };
    const checkMultipleOf = () => {
        if (schemaEntry.multipleOf) {
            // check if default, minimum and maximum are multiples of multipleOf
            if (schemaEntry.default) {
                if (schemaEntry.default % schemaEntry.multipleOf !== 0) {
                    validation['default'] = false;
                } else {
                    validation['default'] = true;
                }
            }
            if (schemaEntry.minimum) {
                if (schemaEntry.minimum % schemaEntry.multipleOf !== 0) {
                    validation['minimum'] = false;
                } else {
                    validation['minimum'] = true;
                }
            }
            if (schemaEntry.maximum) {
                if (schemaEntry.maximum % schemaEntry.multipleOf !== 0) {
                    validation['maximum'] = false;
                } else {
                    validation['maximum'] = true;
                }
            }
        }
    };
</script>

<div class="row border align-items-center pt-2">
    <div class="row align-items-center">
        <div data-dnd-handle class="col-auto border-end">
            <i class="fas fa-grip-vertical"></i>
        </div>

        <div class="input-group col">
            <span class="input-group-text" id="basic-addon1"> <i id="icon" class={schemaEntry.fa_icon}></i></span>
            <div class="form-floating">
                <input type="text" class="form-control" placeholder="icon" value={schemaEntry.fa_icon} />
                <label for="icon">icon</label>
            </div>
        </div>
        <div class="col">
            <slot />
        </div>
        <div class="col">
            <div class="form-floating">
                <select
                    id="type"
                    class="form-select param_key param_type"
                    data-param_key="type"
                    on:change={(e) => updateSchema(e, 'type')}
                >
                    {#each schemaTypes as type}
                        <option selected={type == schemaEntry.type} value={type}>{type}</option>
                    {/each}
                </select>
                <label for="type">type</label>
            </div>
        </div>
        <div class="col-auto d-flex flex-column">
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
    <div class="row">
        <div class="col">
            <div class="form-floating m-3">
                <textarea
                    class="description-textarea form-control"
                    id="description"
                    placeholder="description"
                    data-param_key="description"
                    on:change={(e) => updateSchema(e, 'description')}
                    on:keydown={(e) => updateSchema(e, 'description')}>{sanitizeHtml(schemaEntry.description)}</textarea
                >
                <label for="description">description</label>
            </div>
        </div>
    </div>
    <div class="row">
        <div class="col">
            <div class="form-floating m-3">
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
    </div>
    <div class="d-sm-none w-100"></div>

    <div class="col">
        <div class="form-floating">
            {#if schemaEntry.type == 'boolean'}
                <label for="default">default</label>
                <input
                    type="checkbox"
                    class="form-check-input"
                    class:is-invalid={validation['default'] === false}
                    class:is-valid={validation['default'] === true}
                    id="default"
                    data-param_key="default"
                    value={schemaEntry.default ?? ''}
                    placeholder="default"
                    on:input={(e) => {
                        updateSchema(e, 'default');
                    }}
                />
            {:else}
                <input
                    type={schemaEntry.type == 'integer' || schemaEntry.type == 'number' ? 'number' : 'text'}
                    class="form-control w-auto"
                    class:is-invalid={validation['default'] === false}
                    class:is-valid={validation['default'] === true}
                    id="default"
                    data-param_key="default"
                    value={schemaEntry.default ?? ''}
                    placeholder="default"
                    on:input={(e) => updateSchema(e, 'default')}
                />
                <label for="default">default</label>
            {/if}
        </div>
    </div>
    <div class="col-auto">
        <div class="form-floating">
            <input
                id="enum"
                type="text"
                class="form-control fit-content"
                data-param_key="enum"
                class:is-invalid={validation['enum'] === false}
                class:is-valid={validation['enum'] === true}
                value={schemaEntry.enum ? schemaEntry.enum.join(', ') : ''}
                placeholder="value1, value2, value3"
                on:change={(e) => updateSchema(e, 'enum')}
            />
            <label for="enum">enum</label>
        </div>
    </div>
    {#if schemaEntry.type == 'string'}
        <div class="col-auto">
            <div class="form-floating">
                <input
                    id="pattern"
                    type="text"
                    class="form-control"
                    data-param_key="required"
                    value={schemaEntry.pattern ?? ''}
                    placeholder="pattern"
                    on:change={(e) => {
                        updateSchema(e, 'pattern');
                    }}
                />
                <label for="pattern">pattern</label>
            </div>
            <small class="form-text text-muted">
                Regular expression, used to validate the input string
                <a
                    href="https://json-schema.org/understanding-json-schema/reference/string.html#regular-expressions"
                    target="_blank"
                    data-bs-toggle="tooltip"
                    title="See help about JSON schema regular expressions"
                >
                    <i class="fas fa-question-circle"></i>
                </a>
            </small>
        </div>
        <div class="col-auto">
            <div class="form-floating">
                <select
                    id="format"
                    class="form-select param_key param_format"
                    data-param_key="format"
                    on:change={(e) => updateSchema(e, 'format')}
                >
                    <option value="">none</option>
                    {#each ['date-time', 'time', 'date', 'duration', 'email'] as format}
                        <option selected={format == schemaEntry.format} value={format}>{format}</option>
                    {/each}
                </select>
                <label for="type">format</label>
            </div>
            <small class="form-text text-muted">
                Special format of the string, e.g., date or e-mail, used to validate the input string
                <a
                    href="https://json-schema.org/understanding-json-schema/reference/string#format"
                    target="_blank"
                    data-bs-toggle="tooltip"
                    title="See help about JSON schema formats"
                >
                    <i class="fas fa-question-circle"></i>
                </a>
            </small>
        </div>
    {/if}
    {#if schemaEntry.type == 'number' || schemaEntry.type == 'integer'}
        <div class="col-auto">
            <div class="form-floating">
                <input
                    id="minimum"
                    type="number"
                    class="form-control"
                    class:is-invalid={validation['minimum'] === false}
                    class:is-valid={validation['minimum'] === true}
                    data-param_key="minimum"
                    value={schemaEntry.minimum ?? ''}
                    placeholder="minimum"
                    step={schemaEntry.type == 'integer' ? 1 : schemaEntry.multipleOf ?? 0.1}
                    on:change={(e) => updateSchema(e, 'minimum')}
                />
                <label for="minimum">minimum</label>
            </div>
        </div>
        <div class="col-auto">
            <div class="form-floating">
                <input
                    id="maximum"
                    type="number"
                    class="form-control"
                    data-param_key="maximum"
                    class:is-invalid={validation['maximum'] === false}
                    class:is-valid={validation['maximum'] === true}
                    min={schemaEntry.minimum ?? ''}
                    max={schemaEntry.maximum ?? ''}
                    value={schemaEntry.maximum ?? ''}
                    placeholder="maximum"
                    step={schemaEntry.type == 'integer' ? 1 : schemaEntry.multipleOf ?? 0.1}
                    on:change={(e) => updateSchema(e, 'maximum')}
                />
                <label for="maximum">maximum</label>
            </div>
        </div>
        <div class="col-auto">
            <div class="form-floating">
                <input
                    id="multipleOf"
                    type="number"
                    class="form-control"
                    class:is-invalid={validation['multipleOf'] === false}
                    class:is-valid={validation['multipleOf'] === true}
                    min={schemaEntry.type == 'integer' ? 1 : ''}
                    data-param_key="multipleOf"
                    value={schemaEntry.multipleOf ?? ''}
                    placeholder="multipleOf"
                    on:change={(e) => updateSchema(e, 'multipleOf')}
                />
                <label for="multipleOf">multipleOf</label>
            </div>
        </div>
    {/if}

    <slot name="delete"></slot>
</div>

<style lang="scss">
    .description-textarea,
    .help_text-textarea {
        height: 5rem;
    }
    .form-control:focus {
        width: auto;
    }
</style>
