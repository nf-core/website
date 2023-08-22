<script lang="ts">
    // Import GitHubProfilePicture.svelte component
    import GitHubProfilePicture from './GitHubProfilePicture.svelte';

    export let username: string;
    export let affiliation: string = '';

    const affiliation_str =
        affiliation.length > 0
            ? `<div class="small"><abbr class="small d-inline-block text-truncate" data-bs-placement="bottom" data-bs-toggle="tooltip" title="${affiliation}" style="max-width: 12rem;">${affiliation}</abbr></div>`
            : '';
    const extra_classes = affiliation.length > 0 ? 'small' : '';
</script>

<GitHubProfilePicture
    wrapperClasses="btn btn-light rounded-pill mb-2 p-0 pe-2 d-flex align-items-center"
    image={'https://github.com/' + username + '.png'}
    name={username}
    circle={true}
    size={50}
>
    <div class={'ms-2 d-inline-block text-start d-flex flex-column align-items-center px-2 ' + extra_classes}>
        <slot />
        {@html affiliation_str}
    </div>
</GitHubProfilePicture>
