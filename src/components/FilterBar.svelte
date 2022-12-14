<script>
    import { ButtonGroup, Button, Input } from 'sveltestrap';
    import { onMount } from 'svelte';
    import { createEventDispatcher } from 'svelte';
    import { CurrentFilter, SortBy } from './filter.js';

    export const input = [];
    export let filter;
    export let sortBy;
    export let displayStyle;

    let search = '';
    let dispatch = createEventDispatcher();

    function handleSearch() {
        dispatch('search', search);
    }
    function handleFilter(fil) {
        // remove focus from button
        event.target.blur();

        if ($CurrentFilter.includes(fil.toLowerCase().replace(" ","_"))) {
            CurrentFilter.set($CurrentFilter.filter((f) => f !== fil.toLowerCase().replace(" ","_")));
        } else {
            CurrentFilter.set([...$CurrentFilter, fil.toLowerCase().replace(" ","_")]);
        }
    }
    function handleSort(sor){
        SortBy.set(sor);
    }
    onMount(() => {
        CurrentFilter.set(filter.map((fil) => fil.toLowerCase().replace(" ","_")));
        SortBy.set(sortBy[0]);

        dispatch('search', search);
    });
</script>

<div class="filter-bar">
    <div class="d-none d-md-flex w-100 justify-content-center mb-2">
        <Input class="w-25 me-2" bind:value={search} on:keyup={handleSearch} placeholder="Filter" />
        <div class="ms-3 d-flex align-items-center">
            Filter
            <ButtonGroup class="ms-1 filter-buttons">
                {#each filter as fil}
                    <button
                        type="button"
                        class="btn btn-outline-success text-nowrap"
                        class:active="{$CurrentFilter.includes(fil.toLowerCase().replace(" ","_"))}"
                        on:click={() => handleFilter(fil)}
                        on:mouseout={() => event.target.blur()}
                        on:blur={() => event.target.blur()}
                        >
                        {fil}
                    </button>
                {/each}
            </ButtonGroup>
        </div>
        <div class="ms-3 d-flex align-items-center">
            Sort
            <ButtonGroup class="ms-1 sort-buttons">
                {#each sortBy as sor}
                    <input type="radio" class="btn-check" name="options" id={sor} checked={sor===$SortBy} on:click={()=>handleSort(sor)}/>
                    <label class="btn btn-outline-success text-nowrap" for={sor}>{sor}</label>
                {/each}
            </ButtonGroup>
        </div>
        <div class="ms-3 d-flex align-items-center">
            Display
            <ButtonGroup class="ms-1 display-buttons">
                {#each displayStyle as dis}
                    <Button class="text-nowrap" outline color="success" on:click={() => handleDisplay({ detail: dis })}>{dis}</Button>
                {/each}
            </ButtonGroup>
        </div>
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
                        <Button
                            class="text-nowrap"
                            outline
                            color="success"
                            on:click={() => handleFilter({ detail: fil })}>{fil}</Button
                        >
                    {/each}
                </ButtonGroup>
            </div>
            <div class="ms-3 d-flex align-items-center">
                Sort
                <ButtonGroup class="ms-1 sort-buttons">
                    {#each sortBy as sor}
                        <input type="radio" class="btn-check" name="options" id={sor} autocomplete="off" />
                        <label class="btn btn-outline-success" for={sor}>{sor}</label>
                    {/each}
                </ButtonGroup>
            </div>
            <div class="ms-3 d-flex align-items-center">
                Display
                <ButtonGroup class="ms-1 display-buttons">
                    {#each displayStyle as dis}
                        <Button
                            class="text-nowrap"
                            outline
                            color="success"
                            on:click={() => handleDisplay({ detail: dis })}>{dis}</Button
                        >
                    {/each}
                </ButtonGroup>
            </div>
        </div>
    </div>
</div>

<style>

</style>
