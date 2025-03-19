<script lang="ts">
    import type { CollectionEntry } from "astro:content";

    type EventData = CollectionEntry<"events">["data"];
    type EventLocation = NonNullable<EventData["locations"]>[number];

    interface Props {
        urls: EventLocation[] | EventLocation;
        btnClass?: string;
    }

    let { urls = $bindable(), btnClass = "btn-success" }: Props = $props();
    const urlArray = $derived(Array.isArray(urls) ? urls : [urls]);

    const getIcon = (url: EventLocation) => {
        const link = typeof url.links === "string" ? url.links : (url.links?.[0] ?? "");
        if (link.includes("youtu")) {
            return "fab fa-youtube";
        } else if (link.includes("figshare")) {
            return "fa-solid fa-file-alt";
        } else if (link.includes("zoom.us")) {
            return "fa-solid fa-video";
        } else if (link.includes("gather.town")) {
            return "fa-solid fa-users";
        } else {
            return "fa-solid fa-external-link-alt";
        }
    };

    const getLink = (url: EventLocation) => (typeof url.links === "string" ? url.links : (url.links?.[0] ?? ""));
    const getName = (url: EventLocation) => url.name || url.city || "Join";
</script>

{#if urlArray.length === 1}
    <a class={"btn text-nowrap " + btnClass} href={getLink(urlArray[0])}>
        <i class={getIcon(urlArray[0]) + " me-1"} aria-hidden="true"></i>
        Join now
    </a>
{:else if urlArray.length > 1}
    <div class="dropdown btn-group" role="group">
        <button
            class="btn btn-success me-2 dropdown-toggle text-nowrap"
            type="button"
            data-bs-toggle="dropdown"
            aria-expanded="false"
        >
            Join now
        </button>
        <ul class="dropdown-menu">
            {#each urlArray as url}
                {#if getLink(url)}
                    <li>
                        <a class="dropdown-item" href={getLink(url)}>
                            <i class={getIcon(url) + " me-1"} aria-hidden="true"></i>
                            {getName(url)}
                        </a>
                    </li>
                {/if}
            {/each}
        </ul>
    </div>
{/if}
