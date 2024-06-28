import { test, expect } from '@playwright/test';

test.describe.configure({ mode: 'parallel' });
// @ts-check

test('"getting started" card on docs index page', async ({ page }) => {
    await page.goto('/docs');
    await expect.soft(page).toHaveTitle('Docs');
    await expect.soft(page.locator('.card').locator('.h4').nth(0)).toContainText('Getting started');
});

test('content and order of the left sidebar of the docs page', async ({ page }) => {
    await page.goto('/docs');
    await expect
        .soft(page.locator('.sidebar-left .top-level > li > details > summary .group-header').nth(0))
        .toContainText('Usage');
    await expect
        .soft(page.locator('.sidebar-left .top-level > li > details > summary .group-header').nth(1))
        .toContainText('Contributing');
    await expect
        .soft(page.locator('.sidebar-left .top-level > li > details > summary .group-header').nth(2))
        .toContainText('Tutorials');
    await expect
        .soft(page.locator('.sidebar-left .top-level > li > details > summary .group-header').nth(3))
        .toContainText('Guidelines');
    await expect
        .soft(page.locator('.sidebar-left .top-level > li > details > summary .group-header').nth(4))
        .toContainText('Checklists');
    await expect
        .soft(page.locator('.sidebar-left .top-level > li > details > summary .group-header').nth(5))
        .toContainText('nf-core/tools');
});

test('content of the nf-core tools installation page is rendered correctly', async ({ page }) => {
    await page.goto('/docs/nf-core-tools/installation');
    await expect.soft(page).toHaveTitle('Docs: nf-core/tools: Installation');
    await expect.soft(page.locator('.markdown-content')).toContainText('conda install nf-core');

    // check that left sidebar is expanded to the correct section
    await expect.soft(page.locator('.sidebar-left .top-level > li > details').nth(5)).toHaveAttribute('open');
    // check that left sidebar element is highlighted with a right-hand border
    await expect
        .soft(page.locator('.sidebar-left .top-level > li:nth-child(6) > details ul a').nth(0))
        .toHaveCSS('border-inline-end-color', 'rgb(34, 174, 99)');
});
