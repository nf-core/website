import { test, expect } from '@playwright/test';

test.describe.configure({ mode: 'parallel' });
// @ts-check

test('modules listing page', async ({ page }) => {
    await page.goto('/modules');
    await expect.soft(page).toHaveTitle('nf-core/modules');

    // check that the card with the title "abacas" has the text "included in: viralrecon"
    await expect
        .soft(page.locator('.card', { hasText: 'abacas' }).locator('.card-body'))
        .toContainText('Included in: viralrecon');

    // click on the card with the title "abacas"
    await page.locator('.card a', { hasText: 'abacas' }).click();

    // check that the page title is "modules/abacas"
    await expect.soft(page).toHaveTitle('modules/abacas');
});

test('module page', async ({ page }) => {
    await page.goto('/modules/fastqc');
    await expect.soft(page).toHaveTitle('modules/fastqc');

    // check that the right sidebar contains the text "included in airrflow"
    await expect.soft(page.locator('.sidebar')).toContainText(/\s*included in\s*airrflow/);
});
