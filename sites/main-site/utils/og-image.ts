// Shared fonts and response headers for the satori-rendered PNG endpoints
// (og.png, newsletter header.png).
// Fonts are inlined into the bundle with Vite's `?inline` (a base64 data URI) instead of
// being fetched from Google Fonts on every request. Satori can't read woff2, so these are woff.
import mavenproDataUri from "@assets/fonts/og/maven-pro-700.woff?inline";
import interDataUri from "@assets/fonts/og/inter-400.woff?inline";

const decodeDataUri = (uri: string) => Buffer.from(uri.slice(uri.indexOf(",") + 1), "base64");

export const ogFonts = [
    { name: "mavenpro", data: decodeDataUri(mavenproDataUri), weight: 700 as const, style: "normal" as const },
    { name: "inter", data: decodeDataUri(interDataUri), weight: 400 as const, style: "normal" as const },
];

// Netlify doesn't cache function responses at the edge based on plain Cache-Control, so
// the CDN needs its own header. `durable` shares the cached image across edge nodes;
// Netlify invalidates it on every deploy.
export const ogImageHeaders = {
    "Content-Type": "image/png",
    "Cache-Control": "public, max-age=31536000, immutable",
    "Netlify-CDN-Cache-Control": "public, durable, s-maxage=31536000",
};
