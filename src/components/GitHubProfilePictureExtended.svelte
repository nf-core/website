<script lang="ts">
    // Import GitHubProfilePicture.svelte component
    import GitHubProfilePicture from './GitHubProfilePicture.svelte';

    export let username: string;
    export let affiliation: string = '';
    export let size: number = 50;
    export let wrapperClasses: string = '';
    export let labelClasses: string = '';

    const affiliation_str =
        affiliation.length > 0
            ? `<div class="small"><abbr class="small d-block text-truncate" data-bs-placement="bottom" data-bs-toggle="tooltip" title="${affiliation}" style="max-width: 12rem;">${affiliation}</abbr></div>`
            : '';
    labelClasses = affiliation.length > 0 ? labelClasses + ' small' : labelClasses;
</script>

<GitHubProfilePicture
    {wrapperClasses}
    linkClasses="btn btn-light rounded-pill mb-2 p-0 d-flex align-items-center"
    image={'https://github.com/' + username + '.png'}
    name={username}
    circle={true}
    size={Math.max(size, 25)}
>
    <div class={'ms-2 pe-2 text-start d-flex flex-column profile-name ' + labelClasses}>
        <slot>
            @{username}
        </slot>
        {@html affiliation_str}
    </div>
</GitHubProfilePicture>

<style lang="scss">
    .profile-name :global(p) {
        margin-bottom: 0;
    }
</style>
