<script lang="ts">
    import ListingCard from "@components/ListingCard.svelte";
    import TagSection from "@components/TagSection.svelte";

    interface Props {
        component: {
            name: string;
            path: string;
            type: string;
            meta: {
                description: string;
                name: string;
                keywords?: string[];
                modules?: string[];
                components?: string[];
            };
            pipelines?: {
                name: string;
                version: string;
            }[];
            subworkflows?: string[];
        };
        pipelines?: {
            name: string;
            version: string;
        }[];
        subworkflows?: string[];
    }

    let { component }: Props = $props();
    const href = "/" + component.type + "s/" + component.name + "/";
</script>

<ListingCard>
    {#snippet cardHeader()}
        <div>
            <a class="text-decoration-none" {href}>{@html component.name.replace("_", "_<wbr>")} </a>
            <small class="gh-stats text-small"></small>
        </div>
    {/snippet}

    {#snippet cardBody()}
        <div class="d-flex flex-column justify-content-between h-100">
            <p class="description flex-grow-1 mb-3">{component.meta.description}</p>
            {#if component.meta.keywords}
                <TagSection tags={component.meta.keywords} type="keywords" />
            {/if}
            {#if component.pipelines}
                <TagSection
                    tags={component.pipelines.map((x) => x.name)}
                    type="pipelines"
                    maxShown={3}
                    included
                    inline
                />
            {/if}
            {#if component.type !== "module" && component.meta.components}
                <TagSection tags={component.meta.components} type="modules" maxShown={3} includes inline />
            {/if}
            {#if component.subworkflows}
                <TagSection tags={component.subworkflows} type="subworkflows" maxShown={3} included inline />
            {/if}
            {#if component.meta.modules}
                <TagSection tags={component.meta.modules} type="modules" maxShown={3} includes inline />
            {/if}
        </div>
    {/snippet}
</ListingCard>
