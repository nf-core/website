// @ts-ignore: no types
import { html } from 'satori-html';
import satori, { init as initSatori } from 'satori/wasm';
import sharp from 'sharp';
// @ts-ignore: no types
import initYoga from 'yoga-wasm-web/asm';

const YOGA = initYoga();
initSatori(YOGA);

console.log('rendering og image');
export const get: APIRoute = async ({ params, request }) => {
    const searchParams = new URL(request.url).searchParams;
    const args = Object.fromEntries(searchParams);
    let subtitle = args.subtitle || '';
    subtitle = subtitle.indexOf('. ') > -1 ? subtitle.substring(0, subtitle.indexOf('. ')) : subtitle;
    const text = args.title.length > 25 && subtitle ? args.title : args.title + ' - ' + subtitle;
    const category = args.category || 'nf-core';
    const html_string = `
    <div class="container"
        style="
        height: 100%;
        width: 100%;
        display: flex;
        flexDirection: column;
        alignItems: center;
        justifyContent: space-around;
        fontSize: 32px;
        fontWeight: 600;
        color: #F8F9FA;
        backgroundColor: #212529;
        borderTop: 5pt solid #1a9655;
        backgroundImage: url('https://raw.githubusercontent.com/nf-core/website/new-og-img/public/images/og-img-bg.svg');">


        <div style="display: flex;
        width: 75%;
        height: 100%;
        flexDirection: column;
        alignItems: flex-start;
        justifyContent: space-around;
        paddingRight: 1rem;
        paddingTop: 5rem;
        paddingBottom: 0.5rem;">
            <h1 style="color: #F8F9FA;">${text}</h1>
            <div
                style="fontSize: 28px;
                    borderBottomLeftRadius: 16px;
                    borderBottomRightRadius: 16px;
                    borderTopLeftRadius: 16px;
                    borderTopRightRadius: 16px;
                    padding: 1rem;
                    background: #2c2c2c;
                    border: 5pt solid #757575;">
                ${category}
            </div>
        </div>

    </div>
    <style>
    h1 {
        font-size: 64px;
        font-weight: 500;
    }
    </style>`;

    const imageOptions = { site: request.url, width: 1200, height: 630, debug: false };
    const jsx = html(html_string);
    const buffer = await generateImage(jsx, imageOptions);

    return new Response(buffer, {
        status: 200,
        headers: {
            'Content-Type': 'image/png',
            'Cache-Control': 'max-age=31536000, immutable',
        },
    });
};

type ImageOptions = {
    site: string;
    width: number;
    height: number;
    debug?: boolean;
};

async function generateImage(jsx: any, { width, height, debug }: ImageOptions) {
    const mavenpro = await fetch(
        'https://fonts.gstatic.com/s/mavenpro/v32/7Auup_AqnyWWAxW2Wk3swUz56MS91Eww8cLx1nejpBh8CvRBOA.woff',
    ).then((res) => res.arrayBuffer());
    const svg = await satori(jsx, {
        debug: debug,
        width: width,
        height: height,
        fonts: [
            {
                name: 'Roboto',
                data: mavenpro,
                weight: 700,
                style: 'normal',
            },
        ],
    });

    return await sharp(Buffer.from(svg)).png().toBuffer();
}

export const prerender = false;
