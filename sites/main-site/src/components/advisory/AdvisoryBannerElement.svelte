<script lang="ts">
    import { formatDistanceToNow } from "date-fns";
    import { onMount } from "svelte";
    import type { CollectionEntry } from "astro:content";
    interface Props {
        advisories?: CollectionEntry<"advisories">[];
        advisories_time_category?: string;
        advisories_type_classes?: {};
        advisories_type_icons?: {};
    }

    let {
        advisories = $bindable([]),
        advisories_time_category = "",
        advisories_type_classes = {},
        advisories_type_icons = {},
    }: Props = $props();

    let backgroundIcon = $state("");

    const advisories_duration = (advisories) => {
        advisories.data.advisoriesCountDown = formatDistanceToNow(advisories.data.start);
        if (advisories.data.startDate === advisories.data.endDate) {
            advisories.data.duration =
                new Date(advisories.data.start).toLocaleString("en-US", {
                    year: "numeric",
                    month: "short",
                    day: "numeric",
                    hour: "numeric",
                    minute: "numeric",
                    hour12: false,
                }) +
                "-" +
                new Date(advisories.data.end).toLocaleString("en-US", {
                    hour: "numeric",
                    minute: "numeric",
                    hour12: false,
                });
        } else {
            advisories.data.duration =
                new Date(advisories.data.start).toLocaleString("en-US", {
                    year: "numeric",
                    month: "short",
                    day: "numeric",
                    hour: "numeric",
                    minute: "numeric",
                    hour12: false,
                }) +
                " - " +
                new Date(advisories.data.end).toLocaleString("en-US", {
                    year: "numeric",
                    month: "short",
                    day: "numeric",
                    hour: "numeric",
                    minute: "numeric",
                    hour12: false,
                });
        }
    };

    advisories
        .map((advisories) => {
            if (advisories.data.title.toLowerCase().match("bytesize")) {
                advisories.data.type = "bytesize";
            }
            advisories_duration(advisories);
            return advisories;
        })
        .sort((a, b) => {
            return new Date(a.data.start).getTime() - new Date(b.data.start).getTime();
        });
    if (advisories_time_category === "upcoming") {
        backgroundIcon = "fa-alarm-clock";
        advisories = advisories
            .filter((advisories) => {
                let time_window = 1 * 24 * 60 * 60 * 1000;
                let advisories_start_unix = advisories.data.start.getTime();
                if (advisories.data.announcement?.start !== undefined) {
                    advisories_start_unix = new Date(advisories.data.announcement.start).getTime();
                    time_window = 0;
                }
                const advisories_end_unix = advisories.data.end.getTime();

                // increase time window to a week for advisories longer than 5 hours
                if (
                    advisories_end_unix - advisories_start_unix > 5 * 60 * 60 * 1000 &&
                    advisories.data.announcement?.start === undefined
                ) {
                    time_window = 7 * 24 * 60 * 60 * 1000;
                }
                if (advisories.data.start < new Date() && new Date() < advisories.data.end) {
                    // this is not an upcoming, but an ongoing advisories
                    return false;
                }

                if (advisories_start_unix < new Date().getTime() + time_window && new Date().getTime() < advisories_end_unix) {
                    return true;
                }
            })
            .sort((a, b) => {
                return new Date(a.data.start).getTime() - new Date(b.data.start).getTime();
            });
    } else if (advisories_time_category === "ongoing") {
        backgroundIcon = "fa-broadcast-tower";
        advisories = advisories
            .filter((advisories) => {
                return advisories.data.start < new Date() && new Date() < advisories.data.end;
            })
            .sort((a, b) => {
                return new Date(b.data.start).getTime() - new Date(a.data.start).getTime();
            });
    }

    let heading_title = $derived(
        advisories_time_category.charAt(0).toUpperCase() +
            advisories_time_category.slice(1) +
            " advisories" +
            (advisories.length > 1 ? "s" : ""),
    );
    onMount(() => {
        advisories.map((advisories) => {
            advisories_duration(advisories);
            return advisories;
        });
    });
</script>

{#if advisories.length > 0}
    <div class={advisories_time_category + "-advisories advisories-container border-bottom border-black-subtle"}>
        <div>
            <div class="d-none d-lg-flex">
                <div class="col-lg-4 overflow-hidden ps-3 position-relative d-flex flex-column justify-content-center">
                    <h4 class="display-4 p-2 pb-0 mb-0 flex-grow-1">{heading_title}</h4>
                    <i
                        class={`fad ${backgroundIcon} homepage-header-fa-background mt-5 ms-1 ms-xl-5`}
                        aria-hidden="true"
                    ></i>
                </div>
                <div class="flex-grow-1">
                    {#each advisories as advisories (advisories.id)}
                        <div class="w-100 row align-items-center">
                            <div class="col-8 py-lg-2 text-lg-start">
                                <h5 class="pt-2 pb-0 pb-lg-1">
                                    <a href={"advisories/" + advisories.id + "/"} class="text-success text-decoration-none"
                                        >{advisories.data.title}</a
                                    >
                                    <span class="ms-1 my-auto">
                                        <span class={"badge bg-" + advisories_type_classes[advisories.data.type] + " small"}
                                            ><i class={advisories_type_icons[advisories.data.type] + " me-1"} aria-hidden="true"
                                            ></i>
                                            {advisories.data.type}</span
                                        >
                                    </span>
                                </h5>
                                <p class="lead mb-1">
                                    <a href={"advisories/" + advisories.id + "/"} class="text-body text-decoration-none"
                                        >{@html advisories.data.subtitle}</a
                                    >
                                </p>
                                {#if advisories.data.duration}
                                    <p class="mb-1">
                                        <a
                                            href={"advisories/" + advisories.id + "/"}
                                            class="text-secondary-emphasis text-decoration-none"
                                            >{advisories.data.duration}</a
                                        >
                                    </p>
                                {/if}
                            </div>

                            <div class="col-4 py-lg-2 text-start d-flex flex-column align-items-start">
                                {#if advisories_time_category === "upcoming"}
                                    <div class="text-nowrap ps-1">
                                        <h5>Advisory starts in</h5>
                                        <span class="display-6">
                                            {@html advisories.data.advisoriesCountDown}
                                        </span>
                                    </div>
                                    <div class="btn-group my-2" role="group" aria-label="Advisory details">
                                        <a
                                            href={"advisories/" + advisories.id + "/"}
                                            class="btn btn-outline-success text-nowrap"
                                        >
                                            Advisory Details
                                        </a>
                                        <ExportAdvisoryButton frontmatter={advisories.data} />
                                    </div>
                                {/if}
                                {#if advisories_time_category === "ongoing"}
                                    <div class="">
                                        <div class="btn-group" role="group" aria-label="Advisory details">
                                            <a
                                                href={"advisories/" + advisories.id + "/"}
                                                class="btn btn-outline-success text-nowrap">Advisory Details</a
                                            >
                                            {#if Array.isArray(advisories.data?.locations) && advisories.data.locations.length > 0}
                                                <VideoButton urls={advisories.data.locations} />
                                            {/if}
                                        </div>
                                    </div>
                                {/if}
                            </div>
                        </div>
                        <hr class="mx-4 my-0 py-0" />
                    {/each}
                </div>
            </div>
            <div class="d-lg-none">
                <div class="pt-2 pb-1 mb-2 overflow-hidden mainpage-subheader-heading-header bg-body-tertiary">
                    <h5 class="pt-2 font-weight-light text-center text-sucess">{heading_title}</h5>
                </div>
                {#each advisories as advisories (advisories.id)}
                    <div class="text-center">
                        <h4 class="pt-2 pb-0">
                            <a href={"advisories/" + advisories.id + "/"} class="text-success text-decoration-none"
                                >{advisories.data.title}</a
                            >
                        </h4>
                        <p class="d-sm-none mb-1">
                            <a href={"advisories/" + advisories.id + "/"} class="text-body text-decoration-none"
                                >{@html advisories.data.subtitle}</a
                            ><span class={"badge bg-" + advisories_type_classes[advisories.data.type] + " small ms-3"}
                                ><i class={advisories_type_icons[advisories.data.type] + " me-1"} aria-hidden="true"></i>
                                {advisories.data.type}</span
                            >
                        </p>
                        <div class="small mb-1 mx-3 d-flex flex-column">
                            <a
                                href={"advisories/" + advisories.id + "/"}
                                class="text-secondary-emphasis text-decoration-none mb-2">{advisories.data.duration}</a
                            >
                            <div class="btn-group text-nowrap" role="group" aria-label="Advisory details">
                                <a href={"advisories/" + advisories.id + "/"} class="btn btn-outline-success"> Advisory Details </a>
                                {#if advisories_time_category === "upcoming"}
                                    <ExportAdvisoryButton frontmatter={advisories.data} />
                                {/if}
                                {#if advisories_time_category === "ongoing" && Array.isArray(advisories.data?.locations) && advisories.data.locations.length > 0}
                                    <VideoButton urls={advisories.data.locations} />
                                {/if}
                            </div>
                        </div>
                    </div>
                {/each}
            </div>
        </div>
    </div>
{/if}

<style lang="scss">
    .homepage-header-fa-background {
        position: absolute;
        font-size: 14em;
        opacity: 0.2;
    }
    hr:last-child {
        display: none;
    }
</style>
