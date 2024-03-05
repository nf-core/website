<script lang="ts">
    export let name: string;
    export let count: number = 0;
    export let image: string = 'https://github.com/' + name + '.png';
    export let size: number = 60;
    export let circle: boolean = false;
    export let linkClasses = '';
    export let wrapperClasses = '';
    export let imgClasses = '';

    const tooltip = count > 0 ? `${name} (${count} commits)` : name;
    const avatar_url = image.match(/\?/) ? `${image}&s=${size}` : `${image}?s=${size}`;
</script>

<div class={'github-profile d-flex ' + wrapperClasses}>
    <a
        href="https://github.com/{name}"
        class={'text-decoration-none d-block overflow-scroll ' + linkClasses}
        target="_blank"
        rel="noopener noreferrer"
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
            class={'align-self-end ' + imgClasses}
        />
        <div class="profile-name text-nowrap">
            <slot />
        </div>
    </a>
</div>

<style lang="scss">
    @import '../styles/_variables.scss';
    img {
        min-height: calc(var(--size) * 1px);
        height: 100%;
        background-color: $white;
        min-width: calc(var(--size) * 1px);
    }
    .github-profile {
        container-type: inline-size;
        container-name: github-profile;
        width: 100%;
        // min-width: 15rem;
    }
    .github-profile a {
        width: fit-content;
    }
    @container github-profile (width < 7rem) {
        .profile-name {
            display: none;
        }
    }
</style>
