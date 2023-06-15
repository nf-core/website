<script>
    import icons from 'code-icons';
    import 'code-icons/styles.css';
    import { formatDistanceToNow } from 'date-fns';
    import prettyBytes from 'pretty-bytes';
    import { onMount } from 'svelte';

    onMount(() => {
        console.log('mounted');
        item = {
            name: 'test.txt',
            lastModified: new Date(),
            size: 1000,
        };

        const showPreview = () => {
            console.log('show preview');
        };
    });
    export let item;
    const fileIcon = icons.getClass(item.name);
</script>

<tr class="position-relative">
    <td class="d-flex align-items-center">
        {#if fileIcon}
            <span class={fileIcon + ' file-icon me-2'} />
        {:else}
            <i class="fa-regular fa-file me-2" />
        {/if}

        <span>{item.name}</span>
    </td>
    <td>{formatDistanceToNow(item.lastModified)}</td>
    <td class="text-end">{prettyBytes(item.size)}</td>
</tr>

<style>
    tr {
        cursor: pointer;
    }
</style>
