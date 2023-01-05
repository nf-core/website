<script>
    import { ButtonGroup, Button, Input } from 'sveltestrap';
    import { onMount } from 'svelte';
    import { CurrentFilter, SortBy, DisplayStyle, SearchQuery } from './filter.js';
    export let filter;
    export let sortBy;
    export let displayStyle;

    let search = $SearchQuery;

    function handleSearch(q) {
        SearchQuery.set(q.target.value);
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
        CurrentFilter.set(filter.map((fil) => fil.name));
        SortBy.set(sortBy[0]);
        DisplayStyle.set(displayStyle[0].name);
    });
</script>

<div class="filter-bar">
    <div class="d-none d-md-flex w-100 justify-content-center mb-2">
        <Input class="w-25 me-2" bind:value={search} on:keyup={handleSearch} placeholder="Filter" />
        <div class="ms-3 d-flex align-items-center">
            Filter
            <div class="btn-group ms-1 filter-buttons d-flex" role="group" aria-label="Filter listing">
                {#each filter as fil}
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
        {#if sortBy.length > 1}
            <div class="ms-3 d-flex align-items-center">
                Sort
                <ButtonGroup class="ms-1 sort-buttons">
                    {#each sortBy as sor (sor)}
                        <input
                            type="radio"
                            class="btn-check"
                            name="options"
                            id={sor.replace(' ', '-')}
                            checked={sor === $SortBy}
                            on:click={() => handleSort(sor)}
                        />
                        <label class="btn btn-outline-success text-nowrap" for={sor.replace(' ', '-')}>{sor}</label>
                    {/each}
                </ButtonGroup>
            </div>
        {/if}
        {#if displayStyle.length > 1}
            <div class="ms-3 d-flex align-items-center">
                Display
                <ButtonGroup class="ms-1 display-buttons">
                    {#each displayStyle as dis}
                        <button
                            type="button"
                            class="btn btn-outline-success text-nowrap"
                            on:click={() => handleDisplayStyle(dis.name)}
                            class:active={$DisplayStyle === dis.name}
                            ><i class={dis.icon} />
                        </button>
                    {/each}
                </ButtonGroup>
            </div>
        {/if}
    </div>
    <div class="d-md-none dropdown">
        <button class="btn btn-secondary dropdown-toggle" type="button" data-bs-toggle="dropdown" aria-expanded="false">
            Filter & Sort
        </button>
        <div class="dropdown-menu">
            <Input class="me-2" bind:value={search} on:keyup={handleSearch} placeholder="Filter" />
            <div class="ms-3 d-flex flex-colum flex-md-row align-items-center">
                Filter
                <ButtonGroup class="ms-1 filter-buttons">
                    {#each filter as fil}
                        <button
                            type="button"
                            class="btn btn-outline-success text-nowrap"
                            class:active={$CurrentFilter.includes(fil.name)}
                            on:click={() => handleFilter(fil.name)}
                            on:mouseout={() => event.target.blur()}
                            on:blur={() => event.target.blur()}
                        >
                            {fil.name}
                            <span class="badge bg-secondary ms-1">{fil.count}</span>
                        </button>
                    {/each}
                </ButtonGroup>
            </div>
            <div class="ms-3 d-flex align-items-center">
                Sort
                <ButtonGroup class="ms-1 sort-buttons">
                    {#each sortBy as sor}
                        <input
                            type="radio"
                            class="btn-check"
                            name="options"
                            id={sor}
                            checked={sor === $SortBy}
                            on:click={() => handleSort(sor)}
                        />
                        <label class="btn btn-outline-success text-nowrap" for={sor}>{sor}</label>
                    {/each}
                </ButtonGroup>
            </div>
            <div class="ms-3 d-flex align-items-center">
                Display
                <ButtonGroup class="ms-1 display-buttons">
                    {#each displayStyle as dis}
                        <button
                            type="button"
                            class="btn btn-outline-success text-nowrap"
                            on:click={() => handleDisplayStyle(dis.name)}
                            class:active={$DisplayStyle === dis.name}
                            ><i class={dis.icon} />
                        </button>
                    {/each}
                </ButtonGroup>
            </div>
        </div>
    </div>
</div>

<style>

</style>
