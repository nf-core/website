<script lang="ts">
    import { onMount } from 'svelte';

    // import { dnd } from 'svelte-pragmatic-list';
    import SchemaBuilderListingGroup from './SchemaBuilderListingGroup.svelte';
    import SchemaBuilderToolbar from './SchemaBuilderToolbar.svelte';

    import CopyButton from '@components/CopyButton.svelte';
    let schema = {};

    let defName = 'definitions';
    let schema_path = '';

    let status = 'loading';
    let collapseGroups = false;

    // rename key in the schema

    // post schema to the server
    const postSchema = async (postStatus: string) => {
        const response = await fetch(`http://localhost:8000/process-schema`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                message: 'Schema created successfully',
            },
            body: JSON.stringify({ status: postStatus, schema, schema_path }),
        });
        const data = await response.json();
        if (data.status === 'web_builder_edited') {
            status = 'saved';
        }
    };

    onMount(async () => {
        // read schema_path from URL
        const url = new URL(window.location.href);
        schema_path = url.searchParams.get('schema_path') || '';
        if (!schema_path) {
            console.error('No schema_path found in the URL');
            return;
        }
        try {
            const response = await fetch(`http://localhost:8000/process-schema?schema_path=${schema_path}`);
            const data = await response.json();

            schema = data.data;
            status = 'loaded';
        } catch (error) {
            console.error('Error fetching schema:', error);
            status = 'loaded';
            // Handle the error, maybe show an error message to the user
        }
    });
    const removeGroup = (name) => {
        let result = confirm('Are you sure you want to remove this group?');
        if (!result) return;
        delete schema?.definitions[name];
        schema = { ...schema };
    };
</script>

<div>
    <SchemaBuilderToolbar bind:schema>
        <div class="d-flex">
            <button
                class="btn btn-outline-secondary"
                on:click={() => {
                    collapseGroups = !collapseGroups;
                }}>Collapse Groups</button
            >
        </div>
        <button
            id="finish"
            class="btn btn-primary"
            class:disabled={status !== 'loaded'}
            on:click={() => {
                postSchema('web_builder_edited');
            }}>Finish</button
        >
    </SchemaBuilderToolbar>

    <div class="d-flex flex-column flex-xxl-row">
        <div class="w-75">
            <div class="border rounded-md" class:border-success={status === 'saved'}>
                {#if status === 'loading'}
                    <p>Loading schema...</p>
                {:else if status === 'loaded'}
                    <div class="p-2" class:collapse={collapseGroups}>
                        {#each Object.keys(schema[defName]) as key (key)}
                            <SchemaBuilderListingGroup name={key} {schema}>
                                <button
                                    slot="delete-group"
                                    class="btn btn-outline-danger delete-btn"
                                    on:click={() => removeGroup(key)}
                                >
                                    <i class="fas fa-trash"></i>
                                </button>
                            </SchemaBuilderListingGroup>
                        {/each}
                        <div data-dnd-indicator hidden style="height: 20px; background-color: red; width: 100%;"></div>
                    </div>
                {:else if (status = 'saved')}
                    <p class="m-2">Schema saved! See terminal log for more details.</p>
                {/if}
            </div>
        </div>
        <figure class="w-25">
            <figcaption
                data-rehype-pretty-code-title=""
                data-language="json"
                class="d-flex align-items-center justify-content-between"
            >
                <span><i class="fa-solid fa-database ms-1 me-2"></i>nf-params.json</span><CopyButton
                    text={JSON.stringify(schema, null, 2)}
                    label="Copy Schema <i class='fa-regular fa-clipboard'></i>"
                    copiedLabel="Copied Schema <i class='fa-regular px-1 fa-clipboard-check' />"
                    classes={'m-2  copy-code-button btn btn-sm btn-outline-secondary opacity-50'}
                    copiedClasses={'m-2 copy-code-button btn btn-sm btn-outline-success'}
                />
            </figcaption>
            <pre class="text-secondary p-0">
        <code class="" data-language="json">
            {JSON.stringify(schema, null, 2)}
        </code>
    </pre>
        </figure>
    </div>
</div>

<style lang="scss">
    .delete-btn:not(:hover) {
        color: var(--bs-tertiary-color);
        border-color: var(--bs-tertiary-color);
    }
</style>
