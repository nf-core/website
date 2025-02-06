<!-- @migration-task Error while migrating Svelte code: This migration would change the name of a slot making the component unusable -->
<script lang="ts">
    import { onMount } from 'svelte';
    import { CurrentFilter, Filters, SortBy, DisplayStyle, SearchQuery } from '@components/store';

    let {
        filter = [],
        sortBy = [],
        displayStyle = [],
        filter_name,
    } = $props<{
        filter: { name: string; class?: string }[];
        sortBy: string[];
        displayStyle: { name: string; icon: string }[];
        filter_name?: string;
    }>();

    let search = $state($SearchQuery);

    function handleSearch(q: Event) {
        const target = q.target as HTMLInputElement;
        SearchQuery.set(target.value.trim());
    }

    function handleFilter(fil: string) {
        // remove focus from button
        (event.target as HTMLElement).blur();
        if ($CurrentFilter.find((f) => f.name === fil)) {
            CurrentFilter.set($CurrentFilter.filter((f) => f.name !== fil));
        } else {
            CurrentFilter.set([...$CurrentFilter, { name: fil }]);
        }
    }

    function handleExlusiveFilter(fil: string) {
        // remove focus from button
        (event.target as HTMLElement).blur();
        CurrentFilter.set([{ name: fil }]);
    }

    function handleSort(sor: string) {
        SortBy.set(sor);
    }

    function handleDisplayStyle(style: string) {
        DisplayStyle.set(style);
    }

    onMount(() => {
        $inspect(filter);
        $inspect($CurrentFilter);
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
            class="form-control w-25 me-2 searchbar"
            value={search}
            oninput={handleSearch}
            placeholder="&#xf002; Search..."
        />

        {#if $Filters.length > 0 && $Filters[0].name}
            <div class="d-none d-xl-block ms-3 d-flex align-items-center">
                <div class="btn-group ms-1 filter-buttons d-flex" role="group" aria-label="Filter listing">
                    {#each $Filters as fil}
                        <button
                            type="button"
                            data-bs-toggle="tooltip"
                            data-bs-placement="top"
                            data-bs-delay="500"
                            title={"Double click to only show items from this category"}
                            class={fil.class
                                ? "btn text-nowrap flex-fill btn-outline-" + fil.class
                                : "btn text-nowrap w-100 btn-outline-success"}
                            class:active={$CurrentFilter.find((f) => f.name === fil.name)}
                            onclick={() => handleFilter(fil.name)}
                            ondblclick={() => handleExlusiveFilter(fil.name)}
                            onmouseout={() => event.target.blur()}
                            onblur={() => event.target.blur()}
                        >
                            {#if fil.icon}
                                <i class={fil.icon + ' me-1'}></i>
                            {/if}
                            {fil.name}
                            {#if fil.count >= 0}
                                <span class="badge bg-secondary ms-1">{fil.count}</span>
                            {/if}
                        </button>
                    {/each}
                </div>
            </div>
            <div class="d-xl-none ms-1 ms-md-3 align-items-center">
                <div class="dropdown">
                    <button
                        class="btn btn-outline-success dropdown-toggle text-nowrap"
                        type="button"
                        data-bs-toggle="dropdown"
                        aria-expanded="false"
                    >
                        {#if filter_name}
                            {@render filter_name()}
                        {:else}
                            <span>Status</span>
                        {/if}
                    </button>
                    <ul class="dropdown-menu">
                        {#each $Filters as fil}
                            <li>
                                <div
                                    class="dropdown-item"
                                    title="Filter"
                                    class:active={$CurrentFilter.find((f) => f.name === fil.name)}
                                    onclick={() => handleFilter(fil.name)}
                                    onkeydown={() => handleFilter(fil.name)}
                                    role="button"
                                    tabindex="0"
                                >
                                    {#if fil.icon}
                                        <i class={fil.icon + ' fa-fw me-1'}></i>
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
                        <i class="fa-solid fa-arrow-down-wide-short me-1"></i>
                        <span class="d-none d-md-inline">{$SortBy}</span>
                    </button>
                    <ul class="dropdown-menu">
                        {#each sortBy as sor (sor)}
                            <li>
                                <div
                                    class="dropdown-item"
                                    title="sort"
                                    id={sor.replace(" ", "-")}
                                    class:active={sor === $SortBy}
                                    onclick={() => handleSort(sor)}
                                    onkeydown={() => handleSort(sor)}
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
                            onclick={() => handleDisplayStyle(dis.name)}
                            class:active={$DisplayStyle === dis.name}
                            title={dis.name + " view"}
                            data-bs-toggle="tooltip"
                            aria-label={dis.name + ' view'}
                            ><i class={dis.icon}></i>
                        </button>
                    {/each}
                </div>
            </div>
        {/if}
    </div>
</div>

<style lang="scss">
    .searchbar {
        font-family: "Inter Variable", "Inter override", sans-serif, "Font Awesome 6 Pro";
    }
</style>
