<script lang="ts">
    import { onMount } from 'svelte';

    export let progress: number = 0;
    export let label: string | number = '';
    export let size: number = 50;
    export let href: string = '#';
    export let title: string = '';
    export let isCurrent: boolean = false;

    onMount(() => {
        isCurrent = true;
        progress = 1;
        window.onscroll = () => {
            const winScroll = document.body.scrollTop || document.documentElement.scrollTop;
            const height = document.documentElement.scrollHeight - document.documentElement.clientHeight;
            progress = Math.min((winScroll / height) * 100 + 1, 100);
        };
    });
</script>

<a class="text-decoration-none" {href} data-bs-title={title} data-bs-toggle="tooltip">
    <svg width={size} height={size} viewBox="0 0 100 100">
        <g transform="rotate(-90 50 50)">
            <!-- rotate the circle to start from the top -->
            <circle cx={size} cy={size} r={size - 10} fill="none" stroke-width="5" opacity="0.1" />
            <circle
                cx={size}
                cy={size}
                r={size - 10}
                fill="none"
                stroke-width="10"
                stroke-dasharray="283"
                stroke-dashoffset={283 - (283 * progress) / 100}
            />
        </g>
        <text
            x="50%"
            y="50%"
            class:d-none={progress === 100}
            class:fw-bold={isCurrent}
            dominant-baseline="central"
            text-anchor="middle"
            style="fill: var(--bs-body-color)">{label}</text
        >
        <text
            x="50%"
            y="50%"
            class:d-block={progress === 100}
            class:d-none={progress !== 100}
            class:fw-bold={isCurrent}
            class="success"
            dominant-baseline="central"
            text-anchor="middle"
            style="fill: var(--bs-success)">{label}</text
        >
    </svg>
    {title}
</a>

<style>
    circle {
        stroke: var(--bs-success);
    }
    text {
        font-size: 200%;
        &.success {
            fill: var(--bs-success);
        }
    }
</style>
