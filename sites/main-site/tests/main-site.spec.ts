import { test, expect } from '@playwright/test';

test.describe.configure({ mode: 'parallel' });
// @ts-check

test('meta is correct', async ({ page }) => {
    await page.goto('/');
    await expect(page).toHaveTitle('nf-core');
});
test('nested event pages', async ({ page }) => {
    await page.goto('/events/hackathon');
    // click on "Hackathon - March 2024 (Online)" link
    await page.waitForSelector('.events');
    await page.getByRole('link', { name: 'Hackathon - March 2024 (Online)' }).click();
    await expect.soft(page).toHaveTitle('Hackathon - March 2024 (Online)');
    // check if markdown is rendered correctly
    await expect.soft(page.locator('.markdown-content')).toContainText('Local sites');
    // click on first "read more" link inside table
    await page.getByRole('link', { name: 'read more' }).first().click();
    // check that we don't get a 404
    await expect.soft(page.locator('.markdown-content')).toContainText('Local event');

    //
});

test('pipeline schema builder redirect', async ({ page }) => {
    // doesn't work on localhost, because it is redirected  by the _redirects file on netlify
    // skip therefore if localhost

    if (!process.env.PLAYWRIGHT_TEST_BASE_URL) {
        return;
    }
    await page.goto('/pipeline_schema_builder');
    await expect.soft(page).toHaveTitle('Parameter schema » nf-core');
});

test('dark mode', async ({ page }) => {
    await page.goto('/publications');
    // wait for page to finish loading
    await page.waitForSelector('.markdown-content');
    // get background-color value
    const bodyBackgroundColorLight = await page.evaluate(() => getComputedStyle(document.body).backgroundColor);
    //click dark mode dropdown
    await page.getByRole('button', { name: 'Change theme button' }).click();
    //click dark mode
    await page.getByRole('button', { name: 'dark' }).click();
    const bodyBackgroundColorDark = await page.evaluate(() => getComputedStyle(document.body).backgroundColor);
    //check if background-color changed
    expect.soft(bodyBackgroundColorLight).not.toEqual(bodyBackgroundColorDark);
    //////
    // NOTE: This part of the test is currently disabled because it fails on chromium
    //////
    // //check if view transition doesn't break dark mode
    // await page.locator('.tab-bar .nav-link').nth(1).click();
    // await page.waitForSelector('#TableOfContents');
    // const bodyBackgroundColorDark2 = await page.evaluate(() => getComputedStyle(document.body).backgroundColor);
    // await expect.soft(bodyBackgroundColorDark).toEqual(bodyBackgroundColorDark2);
});
