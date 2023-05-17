<script lang="ts">
    import AutoComplete from 'simple-svelte-autocomplete';
    import { onMount } from 'svelte';

    export let possibleResults: { href: string; name: string }[] = [];
    export let placeholder: string = 'Search';

    let selectedItem;

    // interpret enter as a click on the selected item

    // }
    function goToLink(event: KeyboardEvent) {
        if (selectedItem) {
            // click the selected item
            const element = document.getElementById(selectedItem.name);
            if (element) {
                element.click();
            }
        }
    }

    onMount(() => {
        // focus the input on mount
    });
</script>

<AutoComplete
    dropdownClassName="bg-body"
    {placeholder}
    items={possibleResults}
    bind:selectedItem
    labelFieldName="name"
    onChange={goToLink}
>
    <a slot="item" let:item href={item.href} id={item.name} class="text-decoration-none text-body stretched-link"
        >{item.name}</a
    >
</AutoComplete>
