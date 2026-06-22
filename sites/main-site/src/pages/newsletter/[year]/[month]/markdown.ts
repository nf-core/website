import type { APIRoute } from "astro";
import { getCollection } from "astro:content";
import pipelines_json from "@public/pipelines.json";
import { getNewsletterContentData, getNewsletterStaticPathsData } from "@utils/newsletter";
import { renderNewsletterMarkdown } from "@utils/newsletter-render";

const images = import.meta.glob("/src/assets/**");

export async function getStaticPaths() {
    const pipelines = pipelines_json.remote_workflows;
    const { months, allProposals } = await getNewsletterStaticPathsData(getCollection, pipelines);
    return months.map(({ year, month }) => ({
        params: { year: String(year), month: String(month).padStart(2, "0") },
        props: { year, month, allProposals },
    }));
}

export const GET: APIRoute = async ({ props }) => {
    const { year, month, allProposals } = props as { year: number; month: number; allProposals: any[] };
    const pipelines = pipelines_json.remote_workflows;
    const data = await getNewsletterContentData(getCollection, pipelines, year, month, allProposals, images);
    return new Response(renderNewsletterMarkdown(data), {
        headers: { "Content-Type": "text/markdown; charset=utf-8" },
    });
};
