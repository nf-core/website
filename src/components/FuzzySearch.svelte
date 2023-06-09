<script lang="ts">
    import fuzzysort from 'fuzzysort';
    import { onMount } from 'svelte';

    export let possibleResults: { href: string; name: string }[] = [];
    export let placeholder: string = 'Search';

    let results = [];

    function search(event: KeyboardEvent) {
        const input = event.target as HTMLInputElement;
        const value = input.value;
        if (value) {
            // search for the value
            results = fuzzysort.go(value, possibleResults, {
                key: 'name',
            });
        } else {
            results = [];
        }
    }
</script>

<div>
    <input class="form-control form-control" type="text" {placeholder} aria-label={placeholder} on:keyup={search} />
    <div class="search-results dropdown-menu px-2" class:show={results.length > 0}>
        {#if results}
            {#each results as result (result)}
                <li>
                    <a href={result.obj.href} id={result.obj.name} class="text-decoration-none text-body fs-5">
                        {@html fuzzysort.highlight(result, '<span class="text-success">', '</span>')}
                    </a>
                </li>
            {/each}
        {/if}
    </div>
</div>
