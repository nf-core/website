<script lang="ts">
    import type { CollectionEntry } from "astro:content";
    import { formatDistanceToNow } from "date-fns";
    import {
        formatAdvisoryType,
        formatAdvisoryCategory,
        getAdvisoryTypeIcon,
        getAdvisoryTypeClass,
        getAdvisorySeverityIcon,
        getAdvisorySeverityClass,
        getAdvisoryMetadataItems,
    } from "./advisoryUtils";

    export let frontmatter: CollectionEntry<"advisories">["data"];
    export let slug: string = "";
    export let time_category: string = "";
    export let showDescription: boolean = true;
    export let narrow: boolean = false;

    const severity_class = getAdvisorySeverityClass(frontmatter.severity);
    const metadataItems = getAdvisoryMetadataItems(frontmatter);
</script>

<a href={"/advisories/" + slug + "/"} class="advisory-card-link">
    <div
        class={"card mb-3 rounded-0 rounded-end " + frontmatter.severity}
        style="border-left-color:var(--bs-{severity_class});"
    >
        <div class="card-body">
            <div class="card-title">
                <h4 id={"advisories-" + slug.split("/")[1]} class:h5={narrow}>
                    {frontmatter.title}
                </h4>
            </div>
            <div class="card-text">
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
            </div>
        </div>
        <div class="card-footer p-0" class:p-md-2={!narrow} class:d-none={narrow}>
            <div class="d-flex align-items-center justify-content-between">
                <p class="d-none text-wrap mb-0 ms-2 align-middle" class:d-md-inline-block={!narrow}>
                    Published {formatDistanceToNow(new Date(frontmatter.publishedDate))} ago
                </p>
                <div class="d-flex flex-wrap gap-1 mb-3 justify-content-end">
                    {#each frontmatter.type as type}
                        <span class={`badge bg-${getAdvisoryTypeClass(type)} small`}>
                            <i class={`${getAdvisoryTypeIcon(type)} me-1`} aria-hidden="true"></i>
                            {formatAdvisoryType(type)}
                        </span>
                    {/each}
                    <span class={`badge bg-${getAdvisorySeverityClass(frontmatter.severity)} small`}>
                        <i class={`${getAdvisorySeverityIcon(frontmatter.severity)} me-1`} aria-hidden="true"></i>
                        {formatAdvisoryType(frontmatter.severity)}
                    </span>
                </div>
            </div>
        </div>
    </div>
</a>

<style lang="scss">
    @import "bootstrap/scss/functions";
    @import "bootstrap/scss/mixins";
    @import "bootstrap/scss/variables";

    .advisory-card-link {
        text-decoration: none;
        color: inherit;
        display: block;

        &:hover {
            color: inherit;
        }

        .card {
            transition: transform 0.2s ease-in-out;

            &:hover {
                transform: scale(1.02);
            }
        }
    }

    .card.rounded-0 {
        border-left: 5px solid;
    }

    .narrow .btn:first-child {
        border-left: 0;
    }

    @include media-breakpoint-up(md) {
        .btn-group.float-end:not(.narrow) {
            .btn:first-child {
                border-radius: var(--bs-border-radius) 0 0 var(--bs-border-radius) !important;
            }
            :global(.btn.dropdown-toggle) {
                border-radius: 0 var(--bs-border-radius) var(--bs-border-radius) 0 !important;
            }
        }
    }

    @include media-breakpoint-down(md) {
        .btn-group.float-end {
            width: 100%;
            .btn:first-child {
                border-left: 0;
            }
        }
    }
</style>
