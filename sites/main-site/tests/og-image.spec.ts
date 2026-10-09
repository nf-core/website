import { test, expect } from "@playwright/test";
import type { APIResponse } from "@playwright/test";

// Smoke tests for the satori-rendered PNG endpoints. Both run on demand, so a broken
// font import or a missing satori/sharp dependency in the function bundle shows up here.
async function expectPng(res: APIResponse, width: number, height: number) {
    expect(res.ok(), `status ${res.status()}`).toBeTruthy();
    expect(res.headers()["content-type"]).toBe("image/png");
    expect(res.headers()["cache-status"]).toContain("Netlify Durable");

    const body = await res.body();
    // PNG signature, then the IHDR chunk holds width and height as big-endian uint32s.
    expect(body.subarray(0, 8)).toEqual(Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]));
    expect(body.readUInt32BE(16)).toBe(width);
    expect(body.readUInt32BE(20)).toBe(height);
}

test("og.png renders with the default title", async ({ request }) => {
    await expectPng(await request.get("/og.png"), 1200, 630);
});

test("og.png renders with title, subtitle and category", async ({ request }) => {
    const params = new URLSearchParams({
        title: "Playwright test",
        subtitle: "Rendered with the bundled fonts. This sentence gets cut off.",
        category: "blog",
    });
    await expectPng(await request.get(`/og.png?${params}`), 1200, 630);
});

test("newsletter header.png renders", async ({ request }) => {
    const listing = await (await request.get("/newsletter")).text();
    const month = listing.match(/\/newsletter\/\d{4}\/\d{2}/);
    expect(month, "expected at least one newsletter month link on /newsletter").not.toBeNull();

    await expectPng(await request.get(`${month![0]}/header.png`), 1920, 1080);
});
