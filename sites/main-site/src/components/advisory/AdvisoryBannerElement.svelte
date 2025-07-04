<script lang="ts">
    import type { CollectionEntry } from "astro:content";
    import { formatDistanceToNow } from "date-fns";
    import { formatAdvisoryType, formatAdvisoryCategory } from "./advisoryUtils";

    interface Props {
        advisories?: CollectionEntry<"advisories">[];
        advisories_time_category?: string;
        advisories_type_classes?: Record<string, string>;
        advisories_type_icons?: Record<string, string>;
        advisory: CollectionEntry<"advisories">;
        advisory_classes?: Record<string, string>;
        advisory_icons?: Record<string, string>;
    }

    let {
        advisories = $bindable([]),
        advisories_time_category = "",
        advisory_classes = {},
        advisory_icons = {},
    }: Props = $props();

    let backgroundIcon = $state("");

    // Sort advisories by published date
    advisories = advisories.sort((a, b) => {
        return new Date(b.data.publishedDate).getTime() - new Date(a.data.publishedDate).getTime();
    });

    if (advisories_time_category === "upcoming") {
        backgroundIcon = "fa-alarm-clock";
        advisories = advisories.filter((advisory) => {
            // Show advisories that are less than 7 days old
            const publishedDate = new Date(advisory.data.publishedDate);
            const sevenDaysAgo = new Date();
            sevenDaysAgo.setDate(sevenDaysAgo.getDate() - 7);
            return publishedDate > sevenDaysAgo;
        });
    } else if (advisories_time_category === "ongoing") {
        backgroundIcon = "fa-broadcast-tower";
        advisories = advisories.filter((advisory) => {
            // Show advisories that are older than 7 days
            const publishedDate = new Date(advisory.data.publishedDate);
            const sevenDaysAgo = new Date();
            sevenDaysAgo.setDate(sevenDaysAgo.getDate() - 7);
            return publishedDate <= sevenDaysAgo;
        });
    }

    let heading_title = $derived(
        advisories_time_category.charAt(0).toUpperCase() +
            advisories_time_category.slice(1) +
            " advisories" +
            (advisories.length > 1 ? "s" : ""),
    );
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
                    {#each advisories as advisory (advisory.id)}
                        <div class="w-100 row align-items-center">
                            <div class="col-8 py-lg-2 text-lg-start">
                                <h5 class="pt-2 pb-0 pb-lg-1">
                                    <a
                                        href={"advisories/" + advisory.id + "/"}
                                        class="text-success text-decoration-none">{advisory.data.title}</a
                                    >
                                    <span class="ms-1 my-auto">
                                        {#each advisory.data.type as type}
                                            <span class={`badge bg-${advisory_classes[type]} small me-1`}>
                                                <i class={`${advisory_icons[type]} me-1`} aria-hidden="true"></i>
                                                {formatAdvisoryType(type)}
                                            </span>
                                        {/each}
                                        <span class={`badge bg-${advisory_classes[advisory.data.severity]} small`}>
                                            <i
                                                class={`${advisory_icons[advisory.data.severity]} me-1`}
                                                aria-hidden="true"
                                            ></i>
                                            {formatAdvisoryCategory(advisory.data.severity)}
                                        </span>
                                    </span>
                                </h5>
                                <p class="lead mb-1">
                                    <a href={"advisories/" + advisory.id + "/"} class="text-body text-decoration-none"
                                        >{@html advisory.data.subtitle}</a
                                    >
                                </p>
                                <p class="mb-1">
                                    <a
                                        href={"advisories/" + advisory.id + "/"}
                                        class="text-secondary-emphasis text-decoration-none"
                                        >Published {formatDistanceToNow(new Date(advisory.data.publishedDate))} ago</a
                                    >
                                </p>
                            </div>

                            <div class="col-4 py-lg-2 text-start d-flex flex-column align-items-start">
                                <div class="btn-group my-2" role="group" aria-label="Advisory details">
                                    <a
                                        href={"advisories/" + advisory.id + "/"}
                                        class="btn btn-outline-success text-nowrap"
                                    >
                                        Advisory Details
                                    </a>
                                </div>
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
                {#each advisories as advisory (advisory.id)}
                    <div class="text-center">
                        <h4 class="pt-2 pb-0">
                            <a href={"advisories/" + advisory.id + "/"} class="text-success text-decoration-none"
                                >{advisory.data.title}</a
                            >
                        </h4>
                        <p class="d-sm-none mb-1">
                            <a href={"advisories/" + advisory.id + "/"} class="text-body text-decoration-none"
                                >{@html advisory.data.subtitle}</a
                            >
                            <span class="ms-1">
                                {#each advisory.data.type as type}
                                    <span class={`badge bg-${advisory_classes[type]} small me-1`}>
                                        <i class={`${advisory_icons[type]} me-1`} aria-hidden="true"></i>
                                        {formatAdvisoryType(type)}
                                    </span>
                                {/each}
                                <span class={`badge bg-${advisory_classes[advisory.data.severity]} small`}>
                                    <i class={`${advisory_icons[advisory.data.severity]} me-1`} aria-hidden="true"></i>
                                    {formatAdvisoryCategory(advisory.data.severity)}
                                </span>
                            </span>
                        </p>
                        <div class="small mb-1 mx-3 d-flex flex-column">
                            <a
                                href={"advisories/" + advisory.id + "/"}
                                class="text-secondary-emphasis text-decoration-none mb-2"
                                >Published {formatDistanceToNow(new Date(advisory.data.publishedDate))} ago</a
                            >
                            <div class="btn-group text-nowrap" role="group" aria-label="Advisory details">
                                <a href={"advisories/" + advisory.id + "/"} class="btn btn-outline-success">
                                    Advisory Details
                                </a>
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
