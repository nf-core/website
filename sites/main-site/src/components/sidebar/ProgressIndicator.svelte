<script lang="ts">
    import { onMount } from "svelte";

    import { Confetti } from "svelte-confetti";

    export let progress: number = 0;
    export let label: string | number = "";
    export let size: number = 40;
    export let strokeWidth: number = 10;
    export let href: string = "#";
    export let title: string = "";
    export let isCurrent: boolean = false;
    export let onScroll: boolean = false;
    export let confetti: boolean = false;

    // Calculate the circumference of the circle
    const circumference = 2 * Math.PI * (size - strokeWidth);
    if (onScroll) {
        onMount(() => {
            isCurrent = true;
            progress = 1;
            window.onscroll = () => {
                const winScroll = document.body.scrollTop || document.documentElement.scrollTop;
                const height = document.documentElement.scrollHeight - document.documentElement.clientHeight;
                progress = Math.min((winScroll / height) * 100 + 1, 100);
            };
        });
    }
</script>

<a class="text-decoration-none" {href} data-bs-title={`${progress}% of items are checked`} data-bs-toggle="tooltip">
    <span>
        <svg width={size} height={size} viewBox={`0 0 ${size * 2} ${size * 2}`}>
            <!-- rotate the circle to start from the top -->
            <g transform={`rotate(-90 ${size} ${size})`}>
                <circle
                    cx={size}
                    cy={size}
                    r={size - strokeWidth}
                    fill="none"
                    stroke-width={strokeWidth}
                    stroke={isCurrent ? "var(--bs-primary-border-subtle)" : "var(--bs-secondary)"}
                />
                <circle
                    class="circle-progress"
                    data-progress={progress}
                    class:is-current={isCurrent}
                    data-size={size}
                    data-strokeWidth={strokeWidth}
                    cx={size}
                    cy={size}
                    r={size - strokeWidth}
                    fill="none"
                    stroke="var(--bs-primary)"
                    stroke-width={strokeWidth}
                    stroke-dasharray={circumference}
                    stroke-dashoffset={circumference - (circumference * progress) / 100}
                />
            </g>
            <text
                x="50%"
                y="50%"
                class:d-none={progress === 100}
                dominant-baseline="central"
                text-anchor="middle"
                style="fill: var(--bs-body-color)">{label}</text
            >
            <text
                x="50%"
                y="50%"
                class:d-block={progress === 100}
                class:d-none={progress !== 100}
                class="success"
                dominant-baseline="central"
                text-anchor="middle"
                style="fill: var(--bs-success)">{label}</text
            >
        </svg>
    </span>

    {title}
</a>
{#if confetti && progress === 100}
    <Confetti rounded={true} />
{/if}

<style>
    /* Define a CSS class for the circle */
    .circle-progress:not(.is-current) {
        transition: stroke-dashoffset 0.5s ease-in-out;
    }

    /* Update the stroke-dashoffset based on the progress */
    .circle-progress[data-progress="100"] {
        stroke-dashoffset: 0;
        stroke: var(--bs-success);
    }
    text {
        font-size: 200%;
        &.success {
            fill: var(--bs-success);
        }
    }
</style>
