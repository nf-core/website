<script lang="ts">
    import { onMount } from "svelte";
    import { CurrentFilter, Filters, SortBy, DisplayStyle, SearchQuery } from "@components/store";

    interface FilterItem {
        name: string;
        displayName?: string;
        class?: string;
        icon?: string;
        count?: number;
    }

    let {
        filter = [],
        sortBy = [],
        displayStyle = [],
        filterName = () => "Status",
    } = $props<{
        filter: FilterItem[];
        sortBy: string[];
        displayStyle: { name: string; icon: string }[];
        filterName?: () => string;
    }>();
    // Initialize filters once on mount
    onMount(async () => {
        // Reset stores if no filters provided
        if (filter.length === 0) {
            Filters.set([]);
            CurrentFilter.set([]);
        } else {
            Filters.set(filter);
        }

        // Reset sort if no options provided
        if (sortBy.length === 0) {
            SortBy.set("");
        } else if (!$SortBy) {
            SortBy.set(sortBy[0]);
        }

        // Reset display style if no options provided
        if (displayStyle.length === 0) {
            DisplayStyle.set("grid");
        }

        // Always reset search
        SearchQuery.set("");
    });

    function handleSearch(q: Event) {
        const target = q.target as HTMLInputElement;
        SearchQuery.set(target.value.trim());
    }

    function handleFilter(fil: string, e: Event) {
        e.preventDefault();
        // remove focus from button
        (e.target as HTMLElement).blur();

        const newFilters = $CurrentFilter.some((f) => f.name === fil)
            ? $CurrentFilter.filter((f) => f.name !== fil)
            : [...$CurrentFilter, { name: fil }];

        CurrentFilter.set(newFilters);
    }

    function handleExclusiveFilter(fil: string, e: Event) {
        e.preventDefault();
        // remove focus from button
        (e.target as HTMLElement).blur();
        CurrentFilter.set([{ name: fil }]);
    }

    function handleSort(sor: string) {
        SortBy.set(sor);
    }

    function handleDisplayStyle(style: string) {
        DisplayStyle.set(style);
    }
</script>

<div class="filter-bar mb-2 px-2">
    <div class="d-flex w-100 justify-content-between justify-content-md-center">
        <input
            type="text"
            class="form-control w-25 me-2 searchbar"
            value={$SearchQuery}
            oninput={handleSearch}
            placeholder="Search"
        />

        <div class="d-none d-xl-block ms-3 d-flex align-items-center">
            <div class="btn-group ms-1 filter-buttons d-flex" role="group" aria-label="Filter listing">
                {#each $Filters as fil}
                    <button
                        type="button"
                        data-bs-toggle="tooltip"
                        data-bs-placement="top"
                        data-bs-delay="500"
                        title="Double click to only show items from this category"
                        class={fil.class
                            ? "btn text-nowrap flex-fill btn-outline-" + fil.class
                            : "btn text-nowrap w-100 btn-outline-success"}
                        class:active={$CurrentFilter.some((f) => f.name === fil.name)}
                        onclick={(e) => handleFilter(fil.name, e)}
                        ondblclick={(e) => handleExclusiveFilter(fil.name, e)}
                    >
                        {#if fil.icon}
                            <i class={fil.icon + " me-1"}></i>
                        {/if}
                        {fil.displayName || fil.name}
                        {#if fil.count !== undefined && fil.count >= 0}
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
                    {filterName()}
                </button>
                <ul class="dropdown-menu">
                    {#each $Filters as fil}
                        <li>
                            <div
                                class="dropdown-item"
                                title="Filter"
                                class:active={$CurrentFilter.some((f) => f.name === fil.name)}
                                onclick={(e) => handleFilter(fil.name, e)}
                                onkeydown={(e) => e.key === "Enter" && handleFilter(fil.name, e)}
                                role="button"
                                tabindex="0"
                            >
                                {#if fil.icon}
                                    <i class={fil.icon + " fa-fw me-1"}></i>
                                {/if}
                                {fil.displayName || fil.name}
                                {#if fil.count !== undefined && fil.count >= 0}
                                    <span class="badge bg-secondary ms-1">{fil.count}</span>
                                {/if}
                            </div>
                        </li>
                    {/each}
                </ul>
            </div>
        </div>
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
                                    onclick={(e) => handleSort(sor)}
                                    onkeydown={(e) => handleSort(sor)}
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
                            onclick={(e) => handleDisplayStyle(dis.name)}
                            class:active={$DisplayStyle === dis.name}
                            title={dis.name + " view"}
                            data-bs-toggle="tooltip"
                            aria-label={dis.name + " view"}
                            ><i class={dis.icon}></i>
                        </button>
                    {/each}
                </div>
            </div>
        {/if}
    </div>
</div>
