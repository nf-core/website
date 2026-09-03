import { test, expect } from "@playwright/test";
import type { APIRequestContext } from "@playwright/test";

// End-to-end smoke tests for the newsletter routes. The month is resolved from
// the listing page so these don't hard-code a date that ages out.
async function latestMonthPath(request: APIRequestContext): Promise<string> {
    const res = await request.get("/newsletter");
    const html = await res.text();
    const match = html.match(/\/newsletter\/\d{4}\/\d{2}/);
    expect(match, "expected at least one newsletter month link on /newsletter").not.toBeNull();
    return match![0];
}

test("listing and month page render", async ({ page }) => {
    await page.goto("/newsletter");
    await expect(page).toHaveTitle("nf-core Newsletter");
    // Sign-up form is wired in.
    await expect(page.locator(".nl-signup-form")).toBeVisible();

    // Into the most recent month.
    await page.locator('a[href^="/newsletter/20"]').first().click();
    await expect(page).toHaveTitle(/^Newsletter - \w+ \d{4}$/);
    // Styled content + at least one section heading + the month nav all render.
    await expect(page.locator("#newsletter-content")).toBeVisible();
    await expect(page.locator("#newsletter-content h2").first()).toBeVisible();
    await expect(page.locator("select.newsletter-month-select").first()).toBeVisible();
});

test("alternate formats (markdown / simple / rss) render", async ({ request }) => {
    const month = await latestMonthPath(request);

    // NB: assert on the body, not the content-type — the GET handler sets
    // text/markdown in `astro dev`, but the route is prerendered to a static file
    // that Netlify serves as text/plain, so the header isn't reliable in prod.
    const md = await request.get(`${month}/markdown`);
    expect(md.ok()).toBeTruthy();
    expect(await md.text()).toMatch(/^# nf-core\/newsletter,/);

    const simple = await request.get(`${month}/simple`);
    expect(simple.ok()).toBeTruthy();
    expect(await simple.text()).toContain("nf-core/newsletter,");

    const rss = await request.get("/newsletter/rss.xml");
    expect(rss.ok()).toBeTruthy();
    const rssBody = await rss.text();
    expect(rssBody).toContain("<rss");
    expect(rssBody).toContain("nf-core Newsletter");
});
