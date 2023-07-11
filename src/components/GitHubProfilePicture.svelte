<script lang="ts">
    export let name: string;
    export let count: number = 0;
    export let image: string;
    export let size: number = 60;
    export let circle: boolean = false;
    export let wrapperClasses = '';

    const tooltip = count > 0 ? `${name} (${count} commits)` : name;
    const avatar_url = image.match(/\?/) ? `${image}&s=${size}` : `${image}?s=${size}`;
</script>

<a
    href="https://github.com/{name}"
    class={'text-decoration-none me-1 mb-1 ' + wrapperClasses}
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
    />
    <slot />
</a>

<style lang="scss">
    img {
        max-height: calc(var(--size) * 1px);
    }
</style>
