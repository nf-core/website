import baseConfig from "../../astro.config.base.mjs";
import icon from "astro-icon";
import { defineConfig } from "astro/config";

// https://astro.build/config
export default defineConfig({
  ...baseConfig,
  build: {
    ...baseConfig.build,
    assetsPrefix: import.meta.env.NETLIFY_SITE_URL || "https://nf-core-modules-subworkflows.netlify.app/",
  },
  integrations: [
    ...baseConfig.integrations.filter((i) => i.name !== "astro-icon"),
    icon({
      iconDir: "../main-site/src/icons",
      include: baseConfig.integrations.find((i) => i.name === "astro-icon")?.include,
    }),
  ],
});
