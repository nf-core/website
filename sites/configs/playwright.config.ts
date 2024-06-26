import { defineConfig, devices } from '@playwright/test';

import baseConfig from '../../playwright.config';

/**
 * See https://playwright.dev/docs/test-configuration.
 */
export default defineConfig({
    ...baseConfig,
    testDir: './tests/',
});
