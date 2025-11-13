import baseConfig from "../../astro.config.base.mjs";
import icon from "astro-icon";
import { defineConfig } from "astro/config";

// https://astro.build/config
export default defineConfig({
    ...baseConfig,
    build: {
        ...baseConfig.build,
        assetsPrefix:
            process.env.CONTEXT === "production" ? "https://nf-core-docs.netlify.app/" : process.env.DEPLOY_PRIME_URL,
    },
    integrations: [
        ...baseConfig.integrations.filter((i) => i.name !== "astro-icon"),
        icon({
            iconDir: "../main-site/src/icons",
            include: baseConfig.integrations.find((i) => i.name === "astro-icon")?.include,
        }),
    ],
});
