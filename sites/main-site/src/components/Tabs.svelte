<script lang="ts">
    import { CurrentTab } from "@components/store";
    interface Props {
        names?: string[];
        icons?: string[];
        children?: import("svelte").Snippet;
    }

    let { names = [], icons = [], children }: Props = $props();

    CurrentTab.set(names[0]);
</script>

<ul class="nav nav-tabs nav-justified" role="tablist">
    {#each names as name, index}
        <li class="nav-item" role="presentation">
            <button
                class="nav-link border-bottom-0"
                class:active={name === $CurrentTab}
                id={`#${name}-tab`}
                type="button"
                role="tab"
                aria-controls={`${name}-tab-pane`}
                aria-selected={name === $CurrentTab}
                onclick={() => CurrentTab.set(name)}
                onkeydown={(e) => {
                    if (e.key === "Enter") {
                        CurrentTab.set(name);
                    }
                }}>{@html icons[index]}{name}</button
            >
        </li>
    {/each}
</ul>
<div class="tab-content" id="myTabContent">
    {@render children?.()}
</div>
