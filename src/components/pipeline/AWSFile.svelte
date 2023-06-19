<script>
    import icons from 'code-icons';
    import 'code-icons/styles.css';
    import { formatDistanceToNow } from 'date-fns';
    import prettyBytes from 'pretty-bytes';
    import { onMount } from 'svelte';

    export let item;
    let showPreview;

    onMount(() => {
        const showPreview = () => {
            console.log('showing preview');
            document.getElementById('file-preview').innerHTML =
                '<iframe src="' + item.href + '" width="100%" height="100%"></iframe>';
            document.getElementById('file-preview').classList.remove('d-none');
        };
    });
    // const fileIcon = icons.getClass(item.name);
</script>

<a href={'?file=' + item.name} class="" on:click={() => showPreview()}>
    <div class="row">
        <span class="col">
            <!-- {#if fileIcon}
            <span class={fileIcon + ' file-icon me-2'} />
        {:else} -->
            <i class="fa-regular fa-file me-2" />
            <!-- {/if} -->

            <span>{item.name}</span>
        </span>
        <span class="col">{formatDistanceToNow(item.lastModified)}</span>
        <span class="col-1 text-end">{prettyBytes(item.size)}</span>
    </div>
</a>
