import { test, expect } from '@playwright/test';

test.describe.configure({ mode: 'parallel' });
// @ts-check

test('pipeline redirect works for /$pipeline-name/results', async ({ page }) => {
    await page.goto('/rnaseq/results/');
    await expect.soft(page.getByRole('link', { name: 'Results' })).toHaveClass('nav-link active');
    // check if SSR works correctly for results
    await page.locator('.list-group-item').nth(1).click();
    await expect.soft(page.locator('.file-browser')).toContainText('fastqc');
    // check if file preview works

    // wait for fastqc/ directory to load
    await page.getByRole('link', { name: 'fastqc/' }).click();
    // wait for file  with suffix fastqc.html to load
    await page
        .getByRole('link', { name: /.*_fastqc\.html$/ })
        .nth(0)
        .click();
    // expect iframe with file preview to be present
    const locator = page.frameLocator('#file-preview iframe');

    // check if file preview is loaded and contains "FastQC Report"
    expect(await locator.locator('html').textContent()).toContain('Per base sequence quality');
});
