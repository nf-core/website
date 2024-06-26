import { defineConfig, devices } from '@playwright/test';

import baseConfig from '../../playwright.config';

/**
 * Read environment variables from file.
 * https://github.com/motdotla/dotenv
 */
// require('dotenv').config();

/**
 * See https://playwright.dev/docs/test-configuration.
 */
export default defineConfig({
    ...baseConfig,
    testDir: './tests/',
    /* Run your local dev server before starting the tests */
    webServer: {
        command: 'npm run dev',
        url: 'http://localhost:4321/',
        timeout: 120 * 1000,
        reuseExistingServer: false,
    },
});
