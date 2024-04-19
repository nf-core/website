<script lang="ts">
    import { CurrentTab } from '@components/store';
    export let names: string[] = [];

    CurrentTab.set(names[0]);
    if (!CurrentTab) {
        console.log('CurrentTab', CurrentTab);
    }
</script>

<ul class="nav nav-tabs" role="tablist">
    {#each names as name, index}
        <li class="nav-item" role="presentation">
            <button
                class="nav-link border-bottom-0"
                class:active={name === $CurrentTab}
                id={`#${name}-tab`}
                data-bs-toggle="tab"
                data-bs-target={`.${name}-tab-pane`}
                type="button"
                role="tab"
                aria-controls={`${name}-tab-pane`}
                aria-selected={name === $CurrentTab}>{name}</button
            >
        </li>
    {/each}
</ul>
<div class="tab-content" id="myTabContent">
    <slot />
</div>
