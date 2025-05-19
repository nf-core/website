<script lang="ts">
    import type { CollectionEntry } from "astro:content";
    import { formatDistanceToNow } from "date-fns";

    export let frontmatter: CollectionEntry<"advisories">["data"];
    export let slug: string = "";
    export let type: string = "";
    export let time_category: string = "";
    export let showDescription: boolean = true;
    export let narrow: boolean = false;
    import { advisories_type_classes, advisories_type_icons } from "./advisoryTypes";

    const type_class = advisories_type_classes[type];
</script>

<div class={"card mb-3 rounded-0 rounded-end " + type} style="border-left-color:var(--bs-{type_class});">
    <div class="card-body">
        <div class="card-title">
            <h4 id={"advisories-" + slug.split("/")[1]} class:h5={narrow}>
                <a class="text-center" class:text-decoration-none={narrow} href={/advisories/ + slug + "/"}>
                    {frontmatter.title}
                </a>
            </h4>
        </div>
        <div class="card-text">
            {#if showDescription}
                <p class="mb-0">{@html frontmatter.subtitle}</p>
            {/if}
            <div
                class="d-flex align-items-center mt-2 flex-wrap justify-content-start"
                class:justify-content-md-end={!narrow}
            >
                <div class="d-flex flex-wrap gap-1">
                    {#each frontmatter.type as type}
                        <span class={"badge bg-" + advisories_type_classes[type] + " small"}
                            ><i class={advisories_type_icons[type] + " me-1"} aria-hidden="true"></i>
                            {type}</span
                        >
                    {/each}
                    <span class={"badge bg-" + advisories_type_classes[frontmatter.severity] + " small"}
                        ><i class={advisories_type_icons[frontmatter.severity] + " me-1"} aria-hidden="true"></i>
                        {frontmatter.severity}</span
                    >
                </div>
            </div>
        </div>
    </div>
    <div class="card-footer p-0" class:p-md-2={!narrow} class:d-none={narrow}>
        <div class="d-flex align-items-center justify-content-between">
            <p class="d-none text-wrap mb-0 ms-2 align-middle" class:d-md-inline-block={!narrow}>
                Published {formatDistanceToNow(new Date(frontmatter.publishedDate))} ago
            </p>
            <div
                class="btn-group float-end"
                class:w-100={narrow}
                class:narrow
                role="group"
                aria-label="See details"
            >
                <a
                    href={"/advisories/" + slug + "/"}
                    class="btn btn-outline-success text-nowrap rounded-start-0"
                    class:rounded-0={["future"].includes(time_category)}>See details</a
                >
            </div>
        </div>
    </div>
</div>

<style lang="scss">
    @import "bootstrap/scss/functions";
    @import "bootstrap/scss/mixins";
    @import "bootstrap/scss/variables";

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
