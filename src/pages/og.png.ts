// @ts-ignore: no types
import HeroBackgroundSrc from '/images/background.png';
// @ts-ignore: no types
import { html } from 'satori-html';
import satori, { init as initSatori } from 'satori/wasm';
import sharp from 'sharp';
// @ts-ignore: no types
import initYoga from 'yoga-wasm-web/asm';


const YOGA = initYoga();
initSatori(YOGA);

export const get: APIRoute = async ({ params, request }) => {
    const searchParams = new URL(request.url).searchParams;
    const args = Object.fromEntries(searchParams);
    const html_string = `<div class="container"
    style="height: 100%;
        width: 100%;
        display: flex;
        flex-direction: column;
        align-items: center;
        justify-content: space-between;
        font-size: 32px;
        font-weight: 600;
        color: #F8F9FA;
        background-color: #212529;

   "
   >
   <div style="display:flex; align-items: center; flex-grow:5;">
  <img src="https://raw.githubusercontent.com/nf-core/nf-co.re/master/public_html/assets/img/logo/nf-core-logo-darkbg.png" width="709" height="165"/>

</div>

<div style="display:flex; flex-direction:column; width:100%; height:30%;">
<svg fill="#fff" xmlns="http://www.w3.org/2000/svg"  width="100%" height="25%"><path d="M 0 300 L 400 300 L 200 00 z"
        fill="#1a9655" /></svg>

  <div style="display:flex;
  width: 100%;
  align-items:stretch;
  justify-content:center;
  flex-grow: 8;
  width: 100%;
  height: 25%;
  background-color:#1a9655;">

  <h1>${args.title}</h1>
  </div>

<span style="font-size: 32px;text-align:center; padding-left: 15rem; padding-right:15rem; position:absolute; margin-top:-7.5rem;">
A community effort to collect a curated set of analysis pipelines built using Nextflow.
</span>
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
        'https://fonts.gstatic.com/s/mavenpro/v32/7Auup_AqnyWWAxW2Wk3swUz56MS91Eww8cLx1nejpBh8CvRBOA.woff'
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