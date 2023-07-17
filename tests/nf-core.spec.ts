import { test, expect } from '@playwright/test';

// @ts-check

test('meta is correct', async ({ page }) => {
    await page.goto('http://localhost:3000/');

    await expect(page).toHaveTitle('nf-core');
});

test('pipeline redirect works for /$pipeliname', async ({ page }) => {
    await page.goto('http://localhost:3000/rnaseq');
    await expect(page).toHaveTitle('nf-core/rnaseq');
    // check if markdown is rendered correctly
    await expect(page.locator('.markdown-content')).toContainText('nf-core/rnaseq is a bioinformatics pipeline');
    // check if results redirect works
    await page.goto('http://localhost:3000/rnaseq/results/');

    await expect(page.getByRole('link', { name: 'Results' })).toHaveClass('nav-link active');
    // check if SSR works correctly for results
    await page.locator('.list-group-item').nth(1).click();
    await expect(page.locator('.file-browser')).toContainText('fastqc');
});

test('random pipeline page', async ({ page }) => {
    await page.goto('http://localhost:3000/');
    await page.getByRole('link', { name: 'View Pipelines' }).click();
    await page.locator('.card').locator('a').first().click();
});
