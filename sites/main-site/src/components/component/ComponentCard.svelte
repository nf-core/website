<script lang="ts">
    import ListingCard from "@components/ListingCard.svelte";
    import TagSection from "@components/TagSection.svelte";
    import { filterByKeyword } from "@utils/search";

    interface Props {
        component: {
            name: string;
            type: string;
            meta: {
                description: string;
                description_rendered?: string;
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
    }

    let { component }: Props = $props();
    const href = $derived("/" + component.type + "s/" + component.name + "/");
    const displayName = $derived(component.name.replaceAll("_", "_<wbr>"));
    const pipelineNames = $derived(component.pipelines?.map((p) => p.name));
</script>

<ListingCard>
    {#snippet cardHeader()}
        <div>
            <a class="text-decoration-none d-flex align-items-center" {href}
                >{@html displayName}
                {#if component.meta.deprecated}
                    <small class="badge text-bg-danger text-small ms-auto">deprecated</small>
                {/if}</a
            >
        </div>
    {/snippet}

    {#snippet cardBody()}
        <div class="d-flex flex-column justify-content-between h-100">
            <div class="description flex-grow-1 mb-3">
                {#if component.meta.description_rendered}
                    {@html component.meta.description_rendered}
                {:else}
                    {component.meta.description}
                {/if}
            </div>
            {#if component.meta.keywords}
                <TagSection tags={component.meta.keywords} type="keywords" onTagClick={filterByKeyword} />
            {/if}
            {#if pipelineNames}
                <TagSection tags={pipelineNames} type="pipelines" maxShown={3} included inline />
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
