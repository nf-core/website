<script lang="ts">
    interface Props {
        tags: string[];
        type: "components" | "modules" | "subworkflows" | "pipelines" | "keywords";
        maxShown?: number;
        includes?: boolean;
        included?: boolean;
        inline?: boolean;
    }

    const { tags, maxShown = 0, type, inline = false, includes = false, included = false }: Props = $props();
    let expanded = $state(false);

    const typeSpecificClasses =
        type === "keywords"
            ? "bg-body-tertiary text-success"
            : type === "pipelines"
              ? "border text-body border-warning-subtle bg-warning-subtle"
              : type === "subworkflows"
                ? "border text-body border-success-subtle bg-success-subtle"
                : "border text-body border-info-subtle bg-info-subtle"; // modules and components share styling

    const additionalTagsCount = tags.length - maxShown;

    function getHref(tag: string) {
        if (type === "pipelines") return `/${tag.replace("/", "_")}/`;
        else if (type === "keywords") return null;
        return `/${type}/${tag.replace("/", "_")}/`;
    }
    function onclick() {
        expanded = !expanded;
    }
</script>

<div class="tag-section text-body-secondary" class:topics={type === "keywords"}>
    {#if includes || included}
        <span class="text-small">{includes ? "Includes:" : "Included in:"}</span>
    {/if}

    {#each expanded || !maxShown ? tags : tags.slice(0, maxShown) as tag}
        <a href={getHref(tag)} class="badge fw-normal me-2 text-decoration-none {typeSpecificClasses}">
            {tag}
        </a>
    {/each}

    {#if additionalTagsCount > 0 && maxShown}
        <span
            class="text-small text-nowrap"
            title={`click to show all ${type}`}
            {onclick}
            onkeydown={onclick}
            tabindex="0"
            role="button"
            data-bs-toggle="tooltip"
            data-bs-delay="500"
        >
            {#if expanded}
                <span class="text-small">hide</span>
            {:else}
                <span class="text-small"
                    >and {additionalTagsCount} more {additionalTagsCount === 1 ? type.slice(0, -1) : type}</span
                >
            {/if}
        </span>
    {/if}
</div>
