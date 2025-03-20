<script lang="ts">
    import VideoButton from "@components/VideoButton.svelte";
    import ExportEventButton from "@components/event/ExportEventButton.svelte";
    import LocalDateTime from "@components/event/LocalDateTime.svelte";
    import type { CollectionEntry } from "astro:content";

    export let frontmatter: CollectionEntry<"events">["data"];
    export let slug: string = "";
    export let type: string = "";
    export let time_category: string = "";
    export let showDescription: boolean = true;
    export let narrow: boolean = false;

    const event_type_classes = {
        bytesize: "success",
        hackathon: "primary",
        talk: "info",
        training: "warning",
    };

    const type_class = event_type_classes[type];
    const isSameDay = frontmatter.startDate === frontmatter.endDate;
</script>

<div class={"card mb-3 rounded-0 rounded-end " + type} style="border-left-color:var(--bs-{type_class});">
    <div class="card-body">
        <div class="card-title">
            <h4 id={"event-" + slug.split("/")[1]} class:h5={narrow}>
                <a class="text-center" class:text-decoration-none={narrow} href={/events/ + slug + "/"}>
                    {frontmatter.title}
                </a>
                {#if time_category === "current" && frontmatter.locations}
                    <div class="float-end d-none d-md-inline">
                        <VideoButton urls={frontmatter.locations} btnClass="btn-danger" />
                    </div>
                {/if}
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
                    <i class="fa-regular fa-calendar me-2" aria-hidden="true"></i>
                    <LocalDateTime date={frontmatter.start} />
                    <span>&nbsp;-&nbsp;</span>
                    {#if isSameDay}
                        <span
                            >{new Date(frontmatter.end).toLocaleString("en-US", {
                                hour: "numeric",
                                minute: "numeric",
                                hour12: false,
                            })}</span
                        >
                    {:else}
                        <LocalDateTime date={frontmatter.end} />
                    {/if}
                </p>
            </div>
        </div>
    </div>
    <div class="card-footer p-0" class:p-md-2={!narrow} class:d-none={narrow}>
        <div class="d-flex align-items-center justify-content-between">
            <p class="d-none text-wrap mb-0 ms-2 align-middle" class:d-md-inline-block={!narrow}>
                <LocalDateTime date={frontmatter.start} />
                <span>&nbsp;-&nbsp;</span>
                {#if isSameDay}
                    <span
                        >{new Date(frontmatter.end).toLocaleString("en-US", {
                            hour: "numeric",
                            minute: "numeric",
                            hour12: false,
                        })}</span
                    >
                {:else}
                    <LocalDateTime date={frontmatter.end} />
                {/if}
            </p>
            <div
                class="btn-group float-end"
                class:w-100={narrow}
                class:narrow
                role="group"
                aria-label="See details or export calendar event"
            >
                <a
                    href={"/events/" + slug + "/"}
                    class="btn btn-outline-success text-nowrap rounded-start-0"
                    class:rounded-0={["future"].includes(time_category)}>See details</a
                >
                {#if time_category === "future"}
                    <ExportEventButton {frontmatter} add_class={"btn-outline-success " + " rounded-top-0"} />
                {/if}
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
