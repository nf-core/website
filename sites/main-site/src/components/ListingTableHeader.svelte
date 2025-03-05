<script lang="ts">
    import { CurrentFilter, Filters, SortBy, DisplayStyle, SearchQuery } from "@components/store";

    interface Props {
        name?: string;
        title?: string;
        textEnd?: boolean;
        textCenter?: boolean;
    }

    let { name = "", title = "Sort by " + name.toLowerCase(), textEnd = false, textCenter = false }: Props = $props();
    let sortInverse = false;

    function handleSort(sor) {
        if (sor === $SortBy) {
            sortInverse = !sortInverse;
        } else {
            sortInverse = false;
        }
        SortBy.set(sortInverse ? sor + ";inverse" : sor);
    }
</script>

<th
    class="text-lg-nowrap sortable"
    class:text-end={textEnd}
    class:text-center={textCenter}
    scope="col"
    data-bs-toggle="tooltip"
    data-bs-delay="500"
    {title}
    onclick={() => handleSort(name)}
>
    <i
        class="fa-arrow-up-arrow-down me-2 fa-swap-opacity"
        class:fa-duotone={$SortBy.startsWith(name)}
        class:fa-regular={!$SortBy.startsWith(name)}
        class:text-muted={!$SortBy.startsWith(name)}
        class:fa-swap-opacity={!$SortBy.endsWith(";inverse")}
    ></i>
    {name}
</th>

<style lang="scss">
    .sortable {
        cursor: pointer;
        &:hover {
            background-color: var(--bs-secondary-bg-subtle);
        }
        :global([data-bs-theme="dark"]) &:hover {
            color: var(--bs-white);
            background-color: var(--bs-tertiary-bg);
        }
    }
</style>
