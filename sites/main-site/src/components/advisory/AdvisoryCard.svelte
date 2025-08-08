<script lang="ts">
    import type { CollectionEntry } from "astro:content";
    import { formatDistanceToNow } from "date-fns";
    import nextflowIcon from "../../icons/logos/nextflow.svg?raw";
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
            <div class="d-flex align-items-center justify-content-between">
                <h4 class="mb-0" id={"advisories-" + slug.split("/")[1]}>
                    {frontmatter.title}
                </h4>
                {#if frontmatter.category}
                    <div class="d-flex align-items-center fs-6 fw-normal text-muted">
                        <strong>Category: </strong>
                        <span class="ms-1">
                            {frontmatter.category.map((cat) => formatAdvisoryCategory(cat)).join(", ")}
                        </span>
                    </div>
                {/if}
            </div>
        {/snippet}

        {#snippet cardBody()}
            {#if showDescription}
                <p class="mb-4">{@html frontmatter.subtitle}</p>
            {/if}
            {#if metadataItems.length > 0}
                <div class="mt-2 small">
                    <strong>Affects:</strong>
                    {#each metadataItems as item}
                        <span class=" d-flex align-items-center mb-1">
                            <span class="d-flex align-items-center">
                                {#if item.label === "Nextflow"}
                                    {@html nextflowIcon}
                                {:else}
                                    <i class={`fas ${item.icon}`} aria-hidden="true"></i>
                                {/if}
                                <span class="ms-1">{item.label}:</span></span
                            >
                            {#each item.value.split(",") as line}
                                <code class="text-muted ms-1">{line}</code>
                            {/each}
                        </span>
                    {/each}
                </div>
            {/if}
        {/snippet}
        {#snippet cardFooter()}
            <div class="d-flex align-items-center justify-content-between">
                <p class="text-wrap mb-0 text-secondary text-small">
                    Published {formatDistanceToNow(new Date(frontmatter.publishedDate))} ago
                </p>
                <div class="d-flex align-items-center flex-wrap gap-1">
                    {#each frontmatter.type as type}
                        <span class={`badge text-bg-${advisoryClasses[type]} small`}>
                            <i class={`${advisoryIcons[type]} me-1`} aria-hidden="true"></i>
                            {formatAdvisoryType(type)}
                        </span>
                    {/each}
                </div>
                <div class="d-flex align-items-center">
                    <span class="text-muted me-1">Severity:</span>
                    <span class={`badge text-bg-${advisoryClasses[frontmatter.severity]}`}>
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
            :global(.card) {
                border-color: var(--bs-tertiary-color);
            }
            :global(.card-header) {
                background-color: var(--bs-tertiary-bg);
            }
        }
    }
</style>
