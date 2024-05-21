<script lang="ts">
    import { CurrentFilter, Filters, SortBy, DisplayStyle, SearchQuery } from '@components/store';

    export let name: string = '';
    export let title: string = 'Sort by ' + name.toLowerCase();
    export let textEnd: boolean = false;
    export let textCenter: boolean = false;
    let sortInverse = false;

    function handleSort(sor) {
        if (sor === $SortBy) {
            sortInverse = !sortInverse;
        } else {
            sortInverse = false;
        }
        SortBy.set(sortInverse ? sor + ';inverse' : sor);
    }
</script>

<th
    class="text-nowrap sortable"
    class:text-end={textEnd}
    class:text-center={textCenter}
    scope="col"
    data-bs-toggle="tooltip"
    data-bs-delay="500"
    {title}
    on:click={() => handleSort(name)}
>
    <i
        class="fa-arrow-up-arrow-down me-2 fa-swap-opacity"
        class:fa-duotone={$SortBy.startsWith(name)}
        class:fa-regular={!$SortBy.startsWith(name)}
        class:text-muted={!$SortBy.startsWith(name)}
        class:fa-swap-opacity={!$SortBy.endsWith(';inverse')}
    />
    {name}
</th>

<style lang="scss">
    @import '../styles/_variables.scss';
    .sortable {
        cursor: pointer;
        &:hover {
            background-color: $secondary-bg-subtle;
        }
        :global([data-bs-theme='dark']) &:hover {
            color: $white;
            background-color: $body-tertiary-bg-dark;
        }
    }
    @include media-breakpoint-down(lg) {
        .text-nowrap {
            white-space: normal !important;
        }
    }
</style>
