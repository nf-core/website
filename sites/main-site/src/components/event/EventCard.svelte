<script lang="ts">
    import VideoButton from "@components/VideoButton.svelte";
    import ExportEventButton from "@components/event/ExportEventButton.svelte";
    import LocalDateTime from "@components/event/LocalDateTime.svelte";
    import type { CollectionEntry } from "astro:content";
    import ListingCard from "../ListingCard.svelte";
    import { eventTypeClasses, type EventType } from "@components/event/eventUtils";

    interface Props {
        frontmatter: CollectionEntry<"events">["data"];
        slug?: string;
        type: EventType;
        time_category?: string;
        showDescription?: boolean;
    }

    let { frontmatter, slug = "", type, time_category = "", showDescription = true }: Props = $props();

    const isSameDay = $derived(frontmatter.startDate === frontmatter.endDate);
    const endTime = $derived(
        new Date(frontmatter.end).toLocaleString("en-US", { hour: "numeric", minute: "numeric", hour12: false }),
    );
</script>

{#snippet dateRange()}
    <LocalDateTime date={frontmatter.start} />
    <span>&nbsp;-&nbsp;</span>
    {#if isSameDay}
        <span>{endTime}</span>
    {:else}
        <LocalDateTime date={frontmatter.end} />
    {/if}
{/snippet}

<!-- the layout adapts to the card's own width (container query), not the viewport, so it also fits narrow columns -->
<div class="event-card">
    <ListingCard
        footer={true}
        cardClass="mb-3 rounded-start-0 type"
        cardStyle={`border-left-color:var(--bs-${eventTypeClasses[type]})`}
    >
        {#snippet cardHeader()}
            <span id={"event-" + slug.split("/")[1]} class="event-title h3">
                <a class="text-center" href={"/events/" + slug + "/"}>
                    {frontmatter.title}
                </a>
                {#if time_category === "current" && frontmatter.locations}
                    <div class="float-end d-none d-md-inline">
                        <VideoButton urls={frontmatter.locations} btnClass="btn-danger" />
                    </div>
                {/if}
            </span>
        {/snippet}
        {#snippet cardBody()}
            {#if showDescription}
                <div class="mb-0">{@html frontmatter.renderedSubtitle}</div>
            {/if}
            <p class="event-date-body pe-3 mt-3 ms-1">
                <i class="fa-regular fa-calendar me-2" aria-hidden="true"></i>
                {@render dateRange()}
            </p>
        {/snippet}
        {#snippet cardFooter()}
            <div class="d-flex align-items-center justify-content-between">
                <p class="event-date-footer text-wrap mb-0">
                    {@render dateRange()}
                </p>
                <div class="event-actions btn-group" role="group" aria-label="See details or export calendar event">
                    <a
                        href={"/events/" + slug + "/"}
                        class="btn btn-outline-success text-nowrap rounded-start-0"
                        class:rounded-0={time_category === "future"}>See details</a
                    >
                    {#if time_category === "future"}
                        <ExportEventButton {frontmatter} add_class="btn-outline-success rounded-top-0" />
                    {/if}
                </div>
            </div>
        {/snippet}
    </ListingCard>
</div>

<style lang="scss">
    .event-card {
        container: event-card / inline-size;
    }
    .event-date-body {
        display: none;
    }

    // compact layout: smaller title, date in the body, full-width buttons
    @container event-card (width < 36rem) {
        .event-title {
            font-size: 1.25rem; // Bootstrap's .h5
            a {
                text-decoration: none;
            }
        }
        .event-date-body {
            display: block;
        }
        .event-date-footer {
            display: none;
        }
        .event-actions {
            width: 100%;
        }
    }
</style>
