<script lang="ts">
    export let headings: {
        text: string;
        slug: string;
        depth: number;
        fa_icon?: string;
    }[];

    import { showHidden, showHelp } from './store.js';

    import SidebarToc from './SidebarToc.svelte';
</script>

<SidebarToc {headings}>
    <div class="btn-group btn-sm" role="group" aria-label="Show/Hide hidden parameters and help texts">
        <button type="button" class="btn btn-outline-secondary text-nowrap btn-sm" on:click={() => showHidden.set(!$showHidden)}
        >{($showHidden?"Hide":"Show")} hidden</button
    >
    <button type="button" class="btn btn-outline-secondary text-nowrap btn-sm" on:click={() => showHelp.set(!$showHelp)}
        >{($showHelp?"Hide":"Show")} all help texts</button
    >
    </div>

</SidebarToc>
