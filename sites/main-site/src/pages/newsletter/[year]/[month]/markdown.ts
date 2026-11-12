import type { APIRoute } from "astro";
import { getCollection } from "astro:content";
import pipelines_json from "@public/pipelines.json";
import { getNewsletterStaticPaths } from "@utils/newsletter";
import { getNewsletterContentData } from "@utils/newsletter-content";
import { renderNewsletterMarkdown } from "@utils/newsletter-render";

const images = import.meta.glob("/src/assets/**");

export async function getStaticPaths() {
    return getNewsletterStaticPaths(getCollection, pipelines_json.remote_workflows);
}

export const GET: APIRoute = async ({ props }) => {
    const { year, month, allProposals } = props as { year: number; month: number; allProposals: any[] };
    const pipelines = pipelines_json.remote_workflows;
    const data = await getNewsletterContentData(getCollection, pipelines, year, month, allProposals, images);
    return new Response(renderNewsletterMarkdown(data), {
        headers: { "Content-Type": "text/markdown; charset=utf-8" },
    });
};
