<script lang="ts">
    import { CurrentTab } from '@components/store';
    export let names: string[] = [];

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
                on:click={() => CurrentTab.set(name)}
                on:keydown={(e) => {
                    if (e.key === 'Enter') {
                        CurrentTab.set(name);
                    }
                }}>{name}</button
            >
        </li>
    {/each}
</ul>
<div class="tab-content" id="myTabContent">
    <slot />
</div>
