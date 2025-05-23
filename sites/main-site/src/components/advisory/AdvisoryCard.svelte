<script lang="ts">
    import type { CollectionEntry } from "astro:content";
    import { formatDistanceToNow } from "date-fns";

    export let frontmatter: CollectionEntry<"advisories">["data"];
    export let slug: string = "";
    export let time_category: string = "";
    export let showDescription: boolean = true;
    export let narrow: boolean = false;
    import { advisories_type_classes, advisories_type_icons } from "./advisoryTypes";

    const severity_class = advisories_type_classes[frontmatter.severity];
</script>

<a href={"/advisories/" + slug + "/"} class="advisory-card-link">
<div class={"card mb-3 rounded-0 rounded-end " + frontmatter.severity} style="border-left-color:var(--bs-{severity_class});">
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
            {#if frontmatter.category || frontmatter.nextflowVersions || frontmatter.nextflowExecutors || frontmatter.softwareDependencies}
                <div class="mt-2 text-muted small">
                    {#if frontmatter.category}
                        <span class="me-3">
                            <i class="fas fa-layer-group me-1" aria-hidden="true"></i>
                            {frontmatter.category.map(cat => cat.charAt(0).toUpperCase() + cat.slice(1)).join(', ')}
                        </span>
                    {/if}
                    {#if frontmatter.nextflowVersions}
                        <span class="me-3">
                            <i class="fas fa-code-branch me-1" aria-hidden="true"></i>
                            Nextflow {frontmatter.nextflowVersions.join(', ')}
                        </span>
                    {/if}
                    {#if frontmatter.nextflowExecutors}
                        <span class="me-3">
                            <i class="fas fa-server me-1" aria-hidden="true"></i>
                            {frontmatter.nextflowExecutors.join(', ')}
                        </span>
                    {/if}
                    {#if frontmatter.softwareDependencies}
                        <span>
                            <i class="fas fa-object-group me-1" aria-hidden="true"></i>
                            {#if Array.isArray(frontmatter.softwareDependencies)}
                                {frontmatter.softwareDependencies.map(dep =>
                                    typeof dep === 'string' ? dep : `${dep.name} ${dep.versions ? '(' + dep.versions.join(', ') + ')' : ''}`
                                ).join(', ')}
                            {:else}
                                {frontmatter.softwareDependencies}
                            {/if}
                        </span>
                    {/if}
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
                    <span class={`badge bg-${advisories_type_classes[type]} small`}>
                        <i class={`${advisories_type_icons[type]} me-1`} aria-hidden="true"></i>
                        {type.split('_').map(word => word.charAt(0).toUpperCase() + word.slice(1)).join(' ')}
                    </span>
                {/each}
                <span class={`badge bg-${advisories_type_classes[frontmatter.severity]} small`}>
                    <i class={`${advisories_type_icons[frontmatter.severity]} me-1`} aria-hidden="true"></i>
                    {frontmatter.severity.charAt(0).toUpperCase() + frontmatter.severity.slice(1)}
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
