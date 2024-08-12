<script lang="ts">
    export let id: string = '';
    export let schema = {};

    // post schema to the server
    const postSchema = async () => {
        const response = await fetch(`/.netlify/functions/process_schema?id=${id}`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                message: 'Schema created successfully',
            },
            body: JSON.stringify({ status: 'processed', schema }),
        });
        return response.json();
    };
</script>

<div class="d-flex justify-content-between border rounded p-2 mb-2">
    <div class="d-flex">
        <button class="btn btn-outline-secondary me-2">Add Field</button>
        <button class="btn btn-outline-secondary me-2">Add Group</button>
    </div>
    <div class="d-flex">
        <button class="btn btn-outline-secondary">Collapse Groups</button>
    </div>
    <button id="finish" class="btn btn-primary" on:click={postSchema}>Finish</button>
    <slot />
</div>
