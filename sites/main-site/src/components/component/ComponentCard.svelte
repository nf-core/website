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
                description_rendered?: string;
                name: string;
                keywords?: string[];
                modules?: string[];
                components?: string[];
                deprecated?: boolean;
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
    const href = $derived("/" + component.type + "s/" + component.name + "/");
</script>

<ListingCard>
    {#snippet cardHeader()}
        <div>
            <a class="text-decoration-none d-flex align-items-center" {href}
                >{@html component.name.replace("_", "_<wbr>")}
                {#if component.meta.deprecated}
                    <small class="badge text-bg-danger text-small ms-auto">deprecated</small>
                {/if}</a
            >
        </div>
    {/snippet}

    {#snippet cardBody()}
        <div class="d-flex flex-column justify-content-between h-100">
            {#if component.meta.description_rendered}
                <div class="description flex-grow-1 mb-3">{@html component.meta.description_rendered}</div>
            {:else}
                <p class="description flex-grow-1 mb-3">{component.meta.description}</p>
            {/if}
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
