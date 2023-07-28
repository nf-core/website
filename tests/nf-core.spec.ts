import { test, expect } from '@playwright/test';


test.describe.configure({ mode: 'parallel' });
// @ts-check

test('meta is correct', async ({ page }) => {
    await page.goto('/');

    await expect.soft(page).toHaveTitle('nf-core');
});

test('pipeline redirect works for /$pipeliname', async ({ page }) => {
    await page.goto('/rnaseq');
    await expect.soft(page).toHaveTitle('nf-core/rnaseq');
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
    await page.goto('/');
    await page.getByRole('link', { name: 'View Pipelines' }).click();
    await page.locator('.card').locator('a').first().click();
    await expect.soft(page.locator('.markdown-content')).toContainText('Citations');
});
