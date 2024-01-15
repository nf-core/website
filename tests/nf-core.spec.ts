import { test, expect } from '@playwright/test';


test.describe.configure({ mode: 'parallel' });
// @ts-check

test('meta is correct', async ({ page }) => {
    await page.goto('/');

    await expect.soft(page).toHaveTitle('nf-core');
});

test('pipeline redirect works for /$pipeliname', async ({ page }) => {
    await page.goto('/rnaseq');
    await expect.soft(page).toHaveTitle('rnaseq: Introduction');
    // check if markdown is rendered correctly
    await expect.soft(page.locator('.markdown-content')).toContainText('nf-core/rnaseq is a bioinformatics pipeline');
    // check if results redirect works
    await page.goto('/rnaseq/results/');
    await expect.soft(page.getByRole('link', { name: 'Results' })).toHaveClass('nav-link active');
    // check if SSR works correctly for results
    await page.locator('.list-group-item').nth(1).click();
    await expect.soft(page.locator('.file-browser')).toContainText('fastqc');
});

test('random pipeline page', async ({ page }) => {
    // check if CTA button works and random (=the pipeline with the newest release) pipeline page is loaded
    await page.goto('/');
    await page.getByRole('link', { name: 'View Pipelines' }).click();
    await page.locator('.card').locator('a').first().click();
    await expect.soft(page.locator('.markdown-content')).toContainText('Citations');
});

test('dark mode', async ({ page }) => {
    await page.goto('/sarek');
    // wait for page to finish loading
    await page.waitForSelector('.markdown-content');
    // get background-color value
    const bodyBackgroundColorLight = await page.evaluate(() => getComputedStyle(document.body).backgroundColor);
    //click dark mode dropdown
    await page.getByRole('button', { name: 'Change theme' }).click();
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
test('nested event pages', async ({ page }) => {
    await page.goto('/events/hackathon');
    // click on "Hackathon - March 2024 (Online)" link
    await page.waitForSelector('.events');
    await page.getByRole('link', { name: 'Hackathon - March 2024 (Online)' }).click();
    await expect.soft(page).toHaveTitle('Hackathon - March 2024 (Online)');
    // check if markdown is rendered correctly
    await expect.soft(page.locator('.markdown-content')).toContainText('Local sites');
    // click on first "read more" link inside table
    await page.getByRole('link', { name: 'read more' }).click();
    // check that we don't get a 404
    await expect.soft(page).not.toHaveTitle('404');
    await expect.soft(page.locator('.markdown-content')).toContainText('Local event');

    //
});
