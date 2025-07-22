<script lang="ts">
    import type { CollectionEntry } from "astro:content";
    import { formatDistanceToNow } from "date-fns";
    import {
        formatAdvisoryType,
        formatAdvisoryCategory,
        advisoryClasses,
        advisoryIcons,
        getAdvisoryMetadataItems,
    } from "./advisoryUtils";
    import ListingCard from "../ListingCard.svelte";

    export let frontmatter: CollectionEntry<"advisories">["data"];
    export let slug: string = "";
    export let showDescription: boolean = true;

    const metadataItems = getAdvisoryMetadataItems(frontmatter);
</script>

<a href={"/advisories/" + slug + "/"} class="card-link text-decoration-none">
    <ListingCard footer={true}>
        {#snippet cardHeader()}
            <h4 id={"advisories-" + slug.split("/")[1]}>
                {frontmatter.title}
            </h4>
        {/snippet}

        {#snippet cardBody()}
            {#if showDescription}
                <p class="mb-4">{@html frontmatter.subtitle}</p>
            {/if}
            {#if frontmatter.category || metadataItems.length > 0}
                <div class="mt-2 text-muted small">
                    {#if frontmatter.category}
                        <span class="me-3">
                            <i class="fas fa-tags me-1" aria-hidden="true"></i>
                            {frontmatter.category.map((cat) => formatAdvisoryCategory(cat)).join(", ")}
                        </span>
                    {/if}
                    {#each metadataItems as item}
                        <span class="me-3">
                            <i class={`fas ${item.icon} me-1`} aria-hidden="true"></i>
                            {item.label}
                            {item.value}
                        </span>
                    {/each}
                </div>
            {/if}
        {/snippet}
        {#snippet cardFooter()}
            <div class="d-flex align-items-center justify-content-between">
                <p class="text-wrap mb-0 text-secondary">
                    Published {formatDistanceToNow(new Date(frontmatter.publishedDate))} ago
                </p>
                <div class="d-flex flex-wrap gap-1">
                    {#each frontmatter.type as type}
                        <span class={`badge text-bg-${advisoryClasses[type]} small`}>
                            <i class={`${advisoryIcons[type]} me-1`} aria-hidden="true"></i>
                            {formatAdvisoryType(type)}
                        </span>
                    {/each}
                    <span class={`badge bg-${advisoryClasses[frontmatter.severity]} small`}>
                        <i class={`${advisoryIcons[frontmatter.severity]} me-1`} aria-hidden="true"></i>
                        {formatAdvisoryType(frontmatter.severity)}
                    </span>
                </div>
            </div>
        {/snippet}
    </ListingCard>
</a>

<style lang="scss">
    .card-link {
        &:hover {
            color: inherit;
        }
    }
</style>
