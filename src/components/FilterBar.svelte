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
        if ($CurrentFilter.includes(fil)) {
            CurrentFilter.set($CurrentFilter.filter((f) => f !== fil));
        } else {
            CurrentFilter.set([...$CurrentFilter, fil]);
        }
    }
    function handleExlusiveFilter(fil) {
        // remove focus from button
        event.target.blur();

        CurrentFilter.set([fil]);
    }
    function handleSort(sor) {
        SortBy.set(sor);
    }

    function handleDisplayStyle(style) {
        DisplayStyle.set(style);
    }

    onMount(() => {
        Filters.set($CurrentFilter);
        if (filter.length > 0) {
            CurrentFilter.set(filter.map((fil) => fil.name));
        }
        if (sortBy.length > 0) {
            SortBy.set(sortBy[0]);
        }
    });
</script>

<div class="filter-bar mb-2">
    <div class="d-none d-lg-flex w-100 justify-content-center">
        <input
            type="text"
            class="form-control w-25 me-2"
            bind:value={search}
            on:keyup={handleSearch}
            placeholder="Search"
        />
        {#if $Filters.length > 0 && $Filters[0].name}
            <div class="ms-3 d-flex align-items-center">
                Show:
                <div class="btn-group ms-1 filter-buttons d-flex" role="group" aria-label="Filter listing">
                    {#each $Filters as fil}
                        <button
                            type="button"
                            data-bs-toggle="tooltip"
                            data-bs-placement="top"
                            data-bs-delay="500"
                            title={'Double click to only show items from this category'}
                            class={fil.class
                                ? 'btn text-nowrap flex-fill btn-outline-' + fil.class
                                : 'btn text-nowrap w-100 btn-outline-success'}
                            class:active={$CurrentFilter.includes(fil.name)}
                            on:click={() => handleFilter(fil.name)}
                            on:dblclick={() => handleExlusiveFilter(fil.name)}
                            on:mouseout={() => event.target.blur()}
                            on:blur={() => event.target.blur()}
                        >
                            {#if fil.icon}
                                <i class={fil.icon + ' me-1'} />
                            {/if}
                            {fil.name}
                            {#if fil.count >= 0}
                                <span class="badge bg-secondary ms-1">{fil.count}</span>
                            {/if}
                        </button>
                    {/each}
                </div>
            </div>
        {/if}
        {#if sortBy.length > 1}
            <div class="ms-3 d-none d-xl-flex align-items-center">
                <span class="text-nowrap">Sort by: </span>
                <div class="btn-group ms-1 sort-buttons" role="group" aria-label="Sort buttons">
                    {#each sortBy as sor (sor)}
                        <input
                            type="radio"
                            class="btn-check"
                            name="sort"
                            id={sor.replace(' ', '-')}
                            checked={$SortBy.startsWith(sor)}
                            on:click={() => handleSort(sor)}
                        />
                        <label class="btn btn-outline-success text-nowrap" for={sor.replace(' ', '-')}>{sor}</label>
                    {/each}
                </div>
            </div>
            <div class="ms-3 d-inline d-xl-none align-items-center">
                <div class="dropdown">
                    <button
                        class="btn btn-outline-success dropdown-toggle text-nowrap"
                        type="button"
                        data-bs-toggle="dropdown"
                        aria-expanded="false"
                    >
                        Sorted by: <span class="text-nowrap d-none d-xl-inline">{$SortBy}</span>
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
            <div class="ms-3 d-flex align-items-center">
                Display:
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
    <div class="d-lg-none dropdown-center text-center">
        <button class="btn btn-secondary dropdown-toggle" type="button" data-bs-toggle="dropdown" aria-expanded="false">
            Filter & Sort
        </button>
        <div class="dropdown-menu">
            <input
                type="text"
                class="form-control mt-0 mb-2"
                bind:value={search}
                on:keyup={handleSearch}
                placeholder="Filter"
            />
            {#if $Filters.length > 0}
                <div class="mx-1 mb-2 align-items-center">
                    <span>Show:</span>
                    <div class="btn-group-vertical mx-1 filter-buttons d-flex" role="group" aria-label="Filter listing">
                        {#each $Filters as fil}
                            <button
                                type="button"
                                data-bs-toggle="tooltip"
                                data-bs-placement="top"
                                title={'Double click to only show items from this category'}
                                class={fil.class
                                    ? 'btn text-nowrap flex-fill btn-outline-' + fil.class
                                    : 'btn text-nowrap w-100 btn-outline-success'}
                                class:active={$CurrentFilter.includes(fil.name)}
                                on:click={() => handleFilter(fil.name)}
                                on:dblclick={() => handleExlusiveFilter(fil.name)}
                                on:mouseout={() => event.target.blur()}
                                on:blur={() => event.target.blur()}
                            >
                                {#if fil.icon}
                                    <i class={fil.icon + ' me-1'} />
                                {/if}
                                {fil.name}
                                {#if fil.count >= 0}
                                    <span class="badge bg-secondary ms-1">{fil.count}</span>
                                {/if}
                            </button>
                        {/each}
                    </div>
                </div>
            {/if}
            {#if sortBy.length > 1}
                <div class="mx-2 mb-2 d-flex align-items-center">
                    Sort:
                    <div class="btn-group-vertical ms-1 sort-buttons" role="group" aria-label="Sort buttons">
                        {#each sortBy as sor (sor)}
                            <input
                                type="radio"
                                class="btn-check"
                                name="sort-sm"
                                id={sor.replace(' ', '-') + '-sm'}
                                checked={sor === $SortBy}
                                on:click={() => handleSort(sor)}
                            />
                            <label class="btn btn-outline-success text-nowrap" for={sor.replace(' ', '-') + '-sm'}
                                >{sor}</label
                            >
                        {/each}
                    </div>
                </div>
            {/if}
            {#if displayStyle.length > 1}
                <div class="mx-2 d-flex align-items-center">
                    Display:
                    <div class="btn-group ms-1 display-buttons" role="group" aria-label="Sort buttons">
                        {#each displayStyle as dis}
                            <button
                                type="button"
                                class="btn btn-outline-success text-nowrap"
                                on:click={() => handleDisplayStyle(dis.name)}
                                class:active={$DisplayStyle === dis.name}
                                ><i class={dis.icon} />
                            </button>
                        {/each}
                    </div>
                </div>
            {/if}
        </div>
    </div>
</div>

<style>
</style>
