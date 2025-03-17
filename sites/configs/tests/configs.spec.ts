import { test, expect } from '@playwright/test';

test.describe.configure({ mode: 'parallel' });
// @ts-check

test('config listing page', async ({ page }) => {
    await page.goto('/configs');
    await expect.soft(page).toHaveTitle('nf-core/configs');
    // check if markdown is rendered correctly
    await expect.soft(page.locator('.table')).toContainText('awsbatch');
    // click on awsbatch link
    await page.getByRole('link', { name: 'awsbatch' }).click();
    await expect.soft(page).toHaveTitle('nf-core/configs: awsbatch');
});

test('config page', async ({ page }) => {
    await page.goto('/configs/awsbatch');
    await expect.soft(page).toHaveTitle('nf-core/configs: awsbatch');
    // check if markdown is rendered correctly
    await expect.soft(page.locator('.markdown-content')).toContainText('awsbatch Configuration');

    // check if config file preview is correctly rendered as a code block
    await expect.soft(page.locator('.config-code pre')).toContainText('AWSBATCH Cloud Profile');
});
