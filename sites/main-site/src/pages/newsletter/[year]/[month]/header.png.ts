// @ts-ignore: no types
import { html } from "satori-html";
import satori from "satori";
import sharp from "sharp";
import type { APIRoute } from "astro";
import { getCollection } from "astro:content";
import pipelines_json from "@public/pipelines.json";
import { getNewsletterContentData } from "@utils/newsletter-content";

export const prerender = false;

const TAGLINE = "The monthly nf-core community newsletter";

// Logo is 2300x1000 (ratio 2.3); render at the card content width.
const LOGO_W = 1200;
const LOGO_H = Math.round(LOGO_W / 2.3);

export const GET: APIRoute = async ({ params, request }) => {
    const year = Number(params.year);
    const month = Number(params.month);

    // We only need the month name + preview summary; pass no proposals/images to
    // keep this lightweight (previewText doesn't depend on either).
    const pipelines = pipelines_json.remote_workflows;
    const data = await getNewsletterContentData(getCollection, pipelines, year, month, [], {});
    const date = `${data.monthName} ${data.year}`;
    const summary = data.previewText;

    // Satori doesn't fetch remote <img> itself, so fetch the logo bytes from this
    // deployment's own origin (/images/logo/...) and inline them as a data URI.
    const origin = new URL(request.url).origin;
    const logoRes = await fetch(`${origin}/images/logo/nf-core-newsletter-darkbg.png`);
    const logoUri = `data:image/png;base64,${Buffer.from(await logoRes.arrayBuffer()).toString("base64")}`;

    const html_string = `
    <div style="display:flex;flex-direction:column;width:1920px;height:1080px;box-sizing:border-box;
        background-color:#212528;border-top:13px solid #22ae63;border-bottom:13px solid #22ae63;
        padding:77px 115px;font-family:'inter';">
        <div style="display:flex;width:100%;justify-content:space-between;align-items:flex-start;">
            <div style="font-family:'inter';font-size:45px;color:#aeb3ba;">${TAGLINE}</div>
            <div style="font-family:'mavenpro';font-weight:700;font-size:58px;color:#cfd3d8;">${date}</div>
        </div>
        <div style="display:flex;flex:1;flex-direction:column;align-items:center;justify-content:flex-start;">
            <img src="${logoUri}" style="display:block;width:${LOGO_W}px;height:${LOGO_H}px;margin-top:80px;" />
            <div style="margin-top:60px;font-family:'mavenpro';font-weight:700;font-size:61px;color:#22ae63;
                text-align:center;">${summary}</div>
        </div>
    </div>`;

    const buffer = await generateImage(html(html_string), { width: 1920, height: 1080 });

    return new Response(buffer, {
        status: 200,
        headers: {
            "Content-Type": "image/png",
            "Cache-Control": "max-age=31536000, immutable",
        },
    });
};

async function generateImage(jsx: any, { width, height }: { width: number; height: number }) {
    const mavenpro = await fetch(
        "https://fonts.gstatic.com/s/mavenpro/v32/7Auup_AqnyWWAxW2Wk3swUz56MS91Eww8cLx1nejpBh8CvRBOA.woff",
    ).then((res) => res.arrayBuffer());
    const inter = await fetch(
        `https://fonts.gstatic.com/s/inter/v13/UcCO3FwrK3iLTeHuS_fvQtMwCp50KnMw2boKoduKmMEVuLyfAZ9hjp-Ek-_EeA.woff`,
    ).then((res) => res.arrayBuffer());
    const svg = await satori(jsx, {
        width,
        height,
        fonts: [
            { name: "mavenpro", data: mavenpro, weight: 700, style: "normal" },
            { name: "inter", data: inter, weight: 400, style: "normal" },
        ],
    });
    return await sharp(Buffer.from(svg)).png().toBuffer();
}
