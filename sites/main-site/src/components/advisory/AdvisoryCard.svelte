<script lang="ts">
    import VideoButton from "@components/VideoButton.svelte";
    import type { CollectionEntry } from "astro:content";

    export let frontmatter: CollectionEntry<"advisories">["data"];
    export let slug: string = "";
    export let type: string = "";
    export let time_category: string = "";
    export let showDescription: boolean = true;
    export let narrow: boolean = false;

    const advisories_type_classes = {
        known_regression: "success",
        incompatibility: "warning",
        security: "alert",
        performance: "success",
        data_corruption: "primary",
        scientific_advice: "secondary",
        other: "secondary",
    };

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
                <p class="text-nowrap text-center text-md-start pe-3 mt-2 ms-1" class:d-md-none={!narrow}>
                ToDo
                </p>
            </div>
        </div>
    </div>
    <div class="card-footer p-0" class:p-md-2={!narrow} class:d-none={narrow}>
        <div class="d-flex align-items-center justify-content-between">
            <p class="d-none text-wrap mb-0 ms-2 align-middle" class:d-md-inline-block={!narrow}>
            ToDo II
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
