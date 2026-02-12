<script lang="ts">
    import ListingTableHeader from "@components/ListingTableHeader.svelte";
    import PaginationNav from "@components/PaginationNav.svelte";
    import TagSection from "@components/TagSection.svelte";
    import ComponentCard from "@components/component/ComponentCard.svelte";
    import { SortBy, DisplayStyle, SearchQuery, currentPage } from "@components/store";

    interface Component {
        name: string;
        path: string;
        type: string;
        meta: {
            name: string;
            description: string;
            keywords: string[];
            components?: string[];
            input: {}[];
            output: {}[];
        };
        pipelines?: {
            name: string;
            version: string;
        }[];
    }

    let { components = [] } = $props();

    let displayStyle = $derived($DisplayStyle);
    let searchQuery = $derived($SearchQuery);
    let sortBy = $derived($SortBy);
    let sortInverse = $derived(sortBy.endsWith(";inverse"));
    let currentPageValue = $derived($currentPage);

    const searchComponents = (component: Component) => {
        if (searchQuery === "") {
            return true;
        }
        if (component.meta.name && component.meta.name.toLowerCase().includes(searchQuery.toLowerCase())) {
            return true;
        }
        if (
            component.meta.description &&
            component.meta.description.toLowerCase().includes(searchQuery.toLowerCase())
        ) {
            return true;
        }
        if (
            component.meta.keywords &&
            component.meta.keywords.some((keyword) => keyword.toLowerCase().includes(searchQuery.toLowerCase()))
        ) {
            return true;
        }
        return false;
    };

    const sortComponents = (a: Component, b: Component) => {
        if (sortBy.startsWith("Name")) {
            if (sortInverse) {
                return b.name.localeCompare(a.name);
            } else {
                return a.name.localeCompare(b.name);
            }
        } else if (sortBy.startsWith("# Pipeline integrations")) {
            const aCount = a.pipelines?.length || 0;
            const bCount = b.pipelines?.length || 0;

            if (sortInverse) {
                return aCount - bCount;
            } else {
                return bCount - aCount;
            }
        }
        return 0;
    };

    let filteredComponents = $derived(components.filter(searchComponents).sort(sortComponents));

    let pageSize = $derived(displayStyle === "grid" ? 12 : 25);
    let lastPage = $derived(Math.ceil(filteredComponents.length / pageSize));
    let paginatedItems = $derived(
        filteredComponents.slice((currentPageValue - 1) * pageSize, currentPageValue * pageSize),
    );

    SearchQuery.subscribe(() => {
        $currentPage = 1;
    });
</script>

<div class={`listing px-0 px-lg-2 py-4 ${components.length > 0 ? components[0].type : ""}`}>
    {#if displayStyle === "grid"}
        <div class="grid">
            {#if filteredComponents.length === 0 && searchQuery !== ""}
                <div class="g-col-12 g-col-md-8 g-start-md-3">
                    <div class="alert alert-secondary text-center" role="alert">
                        No components found. Try changing your search query or filters.
                    </div>
                </div>
            {:else}
                {#each paginatedItems as component (component.name)}
                    <div
                        class={[
                            "g-col-12",
                            "g-col-md-6",
                            "g-col-xl-6",
                            "g-col-xxl-4",
                            components[0].type === "module" && "g-col-xxxl-3",
                            components[0].type === "subworkflow" && "g-col-xxxl-4",
                            components[0].type === "module" && "g-col-xxxxl-2",
                            components[0].type === "subworkflow" && "g-col-xxxxl-3",
                        ]}
                    >
                        <ComponentCard {component} />
                    </div>
                {/each}
            {/if}
        </div>
    {:else}
        <table class="table table-responsive">
            <thead>
                <tr>
                    <ListingTableHeader name="Name" />
                    <th scope="col">Description</th>
                    <th class="keywords" scope="col">Keywords</th>
                    {#if components.length > 0 && components[0].type !== "module"}
                        <th class="components" scope="col">Components</th>
                    {/if}
                    <ListingTableHeader
                        name="# Pipeline integrations"
                        title={"Sort by number of pipelines with " + (components.length > 0 ? components[0].type : "")}
                        textEnd={true}
                    />
                </tr>
            </thead>
            <tbody>
                {#if filteredComponents.length === 0 && searchQuery !== ""}
                    <tr>
                        <td colspan="5" class="text-center">
                            <div class="alert alert-secondary" role="alert">
                                No components found. Try changing your search query or filters.
                            </div>
                        </td>
                    </tr>
                {:else}
                    {#each paginatedItems as component (component.name)}
                        <tr>
                            <td class="name p-0">
                                <div class="position-relative p-3">
                                    <a class="stretched-link" href={"/" + component.type + "s/" + component.name + "/"}
                                        >{@html component.name.replace("_", "_<wbr>")}</a
                                    >
                                </div>
                            </td>
                            <td class="text-small">
                                {component.meta.description}
                            </td>
                            <td class="topics">
                                <TagSection tags={component.meta.keywords} type="keywords" />
                            </td>
                            {#if component.type !== "module"}
                                <td class="components">
                                    <TagSection tags={component.meta.components} type="modules" />
                                </td>
                            {/if}
                            <td class="text-end">
                                {component.pipelines ? component.pipelines.length : 0}
                            </td>
                        </tr>
                    {/each}
                {/if}
            </tbody>
        </table>
    {/if}
</div>
{#if lastPage > 0}
    <PaginationNav {lastPage} />
{/if}

<style>
    .name {
        min-width: 15rem;
        word-break: break-word;
    }
    .keywords {
        max-width: 35rem;
    }
</style>
