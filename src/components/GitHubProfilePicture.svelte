<script lang="ts">
    export let name: string;
    export let count: number = 0;
    export let image: string = 'https://github.com/' + name + '.png';
    export let size: number = 60;
    export let circle: boolean = false;
    export let linkClasses: string = '';
    export let wrapperClasses: string = '';
    export let imgClasses: string = '';
    export let containerQuery: boolean = false;

    const tooltip = count > 0 ? `${name} (${count} commits)` : name;
    const avatar_url = image.match(/\?/) ? `${image}&s=${size}` : `${image}?s=${size}`;
</script>

<div class={wrapperClasses} class:github-profile={containerQuery}>
    <a
        href="https://github.com/{name}"
        class={'text-decoration-none d-block ' + linkClasses}
        target="_blank"
        rel="noopener noreferrer"
        style="--size:{size};"
    >
        <img
            src={avatar_url}
            width={size}
            height={size}
            class:rounded-circle={circle}
            data-bs-placement="bottom"
            data-bs-toggle="tooltip"
            title={tooltip}
            alt={`Github user ${name}`}
            style="--size:{size};"
            class={' ' + imgClasses}
        />
        <div class="profile-name text-nowrap">
            <slot />
        </div>
    </a>
</div>

<style lang="scss">
    @import '../styles/_variables.scss';
    img {
        background-color: $white;
    }
    .github-profile {
        container-type: inline-size;
        container-name: github-profile;
    }
    .github-profile a {
        width: fit-content;
        // height: var(--size);
    }
    @container github-profile (width < 7rem) {
        :global(.profile-name) {
            display: none;
        }
    }
    :global(.github-profile:hover) {
        z-index: 1000;
        :global(.profile-name) {
            display: block !important;
            z-index: 1000;
        }
    }
</style>
