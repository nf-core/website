<script lang="ts">
    import { onMount } from "svelte";

    import { Confetti } from "svelte-confetti";

    interface Props {
        progress?: number;
        label?: string | number;
        size?: number;
        strokeWidth?: number;
        href?: string;
        title?: string;
        isCurrent?: boolean;
        onScroll?: boolean;
        confetti?: boolean;
    }

    let {
        progress = 0,
        label = "",
        size = 40,
        strokeWidth = 10,
        href = "#",
        title = "",
        isCurrent = false,
        onScroll = false,
        confetti = false,
    } = $props();

    // Calculate the circumference of the circle
    const circumference = $derived(2 * Math.PI * (size - strokeWidth));

    onMount(() => {
        if (onScroll) {
            isCurrent = true;
            progress = 1;
            window.onscroll = () => {
                const winScroll = document.body.scrollTop || document.documentElement.scrollTop;
                const height = document.documentElement.scrollHeight - document.documentElement.clientHeight;
                progress = Math.min((winScroll / height) * 100 + 1, 100);
            };
        }
    });
</script>

{#if confetti && progress === 100}
    <Confetti rounded={true} />
{/if}
<div class="progress-indicator d-inline-flex align-items-center">
    <svg width={size} height={size} viewBox={`0 0 ${size * 2} ${size * 2}`}>
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
    {#if title}
        <span class="ms-1">{title}</span>
    {/if}
</div>

<style>
    .circle-progress:not(.is-current) {
        transition: stroke-dashoffset 0.5s ease-in-out;
    }

    .circle-progress[data-progress="100"] {
        stroke-dashoffset: 0;
        stroke: var(--bs-success);
    }

    text {
        font-size: 200%;
    }
</style>
