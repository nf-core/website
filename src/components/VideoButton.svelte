<script>
    export let urls = [];
    export let btnClass = 'btn-success';
    if (typeof urls === 'string') {
        urls = [urls];
    }
    // check if url is for youtube, figshare, zoom, or gather.town and return the appropriate icon
    const getIcon = (url) => {
        if (url.includes('youtu')) {
            return 'fab fa-youtube';
        } else if (url.includes('figshare')) {
            return 'fa-solid fa-file-alt';
        } else if (url.includes('zoom.us')) {
            return 'fa-solid fa-video';
        } else if (url.includes('gather.town')) {
            return 'fa-solid fa-users';
        } else {
            return 'fa-solid fa-external-link-alt';
        }
    };
</script>

{#if urls.length === 1}
    <a class={'btn ' + btnClass} href={urls[0]}>
        <i class={getIcon(urls[0]) + ' me-1'} aria-hidden="true" />
        Join now
    </a>
{:else}
    <div class="dropdown btn-group" role="group">
        <button
            class="btn btn-success me-2 dropdown-toggle"
            type="button"
            id="dropdownMenuButton1"
            data-bs-toggle="dropdown"
            aria-expanded="false"
        >
            Join now
        </button>
        <ul class="dropdown-menu" aria-labelledby="dropdownMenuButton1">
            {#each urls as url (url)}
                <li>
                    <a class="dropdown-item" href={url}>
                        <i class={getIcon(url) + ' me-1'} aria-hidden="true" />
                        {url}
                    </a>
                </li>
            {/each}
        </ul>
    </div>
{/if}
