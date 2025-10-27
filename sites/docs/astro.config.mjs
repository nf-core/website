import baseConfig from "../../astro.config.base.mjs";
import { defineConfig } from "astro/config";

// https://astro.build/config
export default defineConfig({
  ...baseConfig,
  build: {
    ...baseConfig.build,
    assetsPrefix: "https://nf-core-docs.netlify.app/",
  },
  integrations: baseConfig.integrations.map((integration) => {
    // Override icon integration to set iconDir
    if (integration.name === "astro-icon") {
      return integration.hooks.then((hooks) => ({
        ...hooks,
        config: { ...hooks.config, iconDir: "../main-site/src/icons" },
      }));
    }
    return integration;
  }),
});
