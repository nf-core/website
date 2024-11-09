import { test, expect } from '@playwright/test';

test.describe.configure({ mode: 'parallel' });
// @ts-check

test('pipeline redirect works for /$pipeline-name', async ({ page }) => {
    await page.goto('/rnaseq');
    await expect.soft(page).toHaveTitle('rnaseq: Introduction');
    // check if markdown is rendered correctly
    await expect.soft(page.locator('.markdown-content')).toContainText('nf-core/rnaseq is a bioinformatics pipeline');
});

test('random pipeline page', async ({ page }) => {
    // check if CTA button works and random (=the pipeline with the newest release) pipeline page is loaded
    await page.goto('/pipelines');
    await page.locator('.card').locator('a').first().click();
    await expect.soft(page.locator('.markdown-content')).toContainText('Citations');
});
