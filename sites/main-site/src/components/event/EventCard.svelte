<script lang="ts">
    import VideoButton from "@components/VideoButton.svelte";
    import ExportEventButton from "@components/event/ExportEventButton.svelte";
    import LocalDateTime from "@components/event/LocalDateTime.svelte";
    import type { CollectionEntry } from "astro:content";
    import ListingCard from "../ListingCard.svelte";

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

    const isSameDay = frontmatter.startDate === frontmatter.endDate;
</script>

<ListingCard
    footer={true}
    cardClass={`mb-3 rounded-start-0 type`}
    cardStyle={`border-left-color:var(--bs-${event_type_classes[type]})`}
>
    {#snippet cardHeader()}
        <span id={"event-" + slug.split("/")[1]} class={narrow ? "h5" : "h3"}>
            <a class="text-center" class:text-decoration-none={narrow} href={"/events/" + slug + "/"}>
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
    {/snippet}
    {#snippet cardFooter()}
        <div class="d-flex align-items-center justify-content-between">
            <p class="d-none text-wrap mb-0" class:d-md-inline-block={!narrow}>
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
    {/snippet}
</ListingCard>

<style lang="scss">
    @media (max-width: 768px) {
        .text-nowrap {
            white-space: normal !important;
        }
    }
</style>
