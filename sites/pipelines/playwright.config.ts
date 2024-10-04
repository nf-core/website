import { defineConfig } from '@playwright/test';

import baseConfig from '../../playwright.config';

/**
 * See https://playwright.dev/docs/test-configuration.
 */
export default defineConfig({
    ...baseConfig,
    testDir: './tests/',
    /* Run your local dev server before starting the tests */
    // webServer: {
    //     command: 'npm run dev',
    //     url: 'http://localhost:4321/',
    //     timeout: 120 * 1000,
    //     reuseExistingServer: false,
    // },
});
