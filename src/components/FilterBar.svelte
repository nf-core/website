<script lang="ts">
    import { CurrentFilter, Filters, SortBy, DisplayStyle, SearchQuery } from '@components/store';
    import { onMount } from 'svelte';

    export let filter: { name: string; class?: string }[] = [];
    export let sortBy: string[] = [];
    export let displayStyle: { name: string; icon: string }[] = [];

    let search = $SearchQuery;

    function handleSearch(q) {
        SearchQuery.set(q.target.value.trim());
    }
    function handleFilter(fil) {
        // remove focus from button
        event.target.blur();
        if ($CurrentFilter.find((f) => f.name === fil)) {
            CurrentFilter.set($CurrentFilter.filter((f) => f.name !== fil));
        } else {
            CurrentFilter.set([...$CurrentFilter, { name: fil }]);
        }
    }
    function handleSort(sor) {
        SortBy.set(sor);
    }

    function handleDisplayStyle(style) {
        DisplayStyle.set(style);
    }
    onMount(() => {
        if (filter.length > 0 && !$CurrentFilter.length) {
            CurrentFilter.set(filter);
        }
        if (filter.length > 0) {
            Filters.set(filter);
        }
        if (sortBy.length > 0) {
            SortBy.set(sortBy[0]);
        }
    });
</script>

<div class="filter-bar mb-2 px-2">
    <div class="d-flex w-100 justify-content-between justify-content-md-center">
        <input
            type="text"
            class="form-control w-25 me-2"
            bind:value={search}
            on:keyup={handleSearch}
            placeholder="Search"
        />
        {#if $Filters.length > 0 && $Filters[0].name}
            <div class="ms-1 ms-md-3 align-items-center">
                <div class="dropdown">
                    <button
                        class="btn btn-outline-success dropdown-toggle text-nowrap"
                        type="button"
                        data-bs-toggle="dropdown"
                        aria-expanded="false"
                    >
                        <slot name="filter-name">Status</slot>
                    </button>
                    <ul class="dropdown-menu">
                        {#each $Filters as fil}
                            <li>
                                <div
                                    class="dropdown-item"
                                    title="Filter"
                                    class:active={$CurrentFilter.find((f) => f.name === fil.name)}
                                    on:click={() => handleFilter(fil.name)}
                                    on:keydown={() => handleFilter(fil.name)}
                                    role="button"
                                    tabindex="0"
                                >
                                    {#if fil.icon}
                                        <i class={fil.icon + ' fa-fw me-1'} />
                                    {/if}
                                    {fil.name}
                                    {#if fil.count >= 0}
                                        <span class="badge bg-secondary ms-1">{fil.count}</span>
                                    {/if}
                                </div>
                            </li>
                        {/each}
                    </ul>
                </div>
            </div>
        {/if}
        {#if sortBy.length > 1}
            <div class="ms-1 ms-md-3 align-items-center">
                <div class="dropdown">
                    <button
                        class="btn btn-outline-success dropdown-toggle text-nowrap"
                        type="button"
                        data-bs-toggle="dropdown"
                        aria-expanded="false"
                    >
                        <i class="fa-solid fa-arrow-down-wide-short me-1 d-xxl-none"></i>
                        <span class="d-none d-xxl-inline">Sort: </span>
                        <span class="d-none d-xl-inline">{$SortBy}</span>
                    </button>
                    <ul class="dropdown-menu">
                        {#each sortBy as sor (sor)}
                            <li>
                                <div
                                    class="dropdown-item"
                                    title="sort"
                                    id={sor.replace(' ', '-')}
                                    class:active={sor === $SortBy}
                                    on:click={() => handleSort(sor)}
                                    on:keydown={() => handleSort(sor)}
                                    role="button"
                                    tabindex="0"
                                >
                                    {sor}
                                </div>
                            </li>
                        {/each}
                    </ul>
                </div>
            </div>
        {/if}
        {#if displayStyle.length > 1}
            <div class="ms-1 ms-md-3 d-flex align-items-center">
                <div class="btn-group ms-1 display-buttons" role="group" aria-label="Display style">
                    {#each displayStyle as dis}
                        <button
                            type="button"
                            class="btn btn-outline-success text-nowrap"
                            on:click={() => handleDisplayStyle(dis.name)}
                            class:active={$DisplayStyle === dis.name}
                            title={dis.name + ' view'}
                            data-bs-toggle="tooltip"
                            ><i class={dis.icon} />
                        </button>
                    {/each}
                </div>
            </div>
        {/if}
    </div>
</div>

<style>
</style>
