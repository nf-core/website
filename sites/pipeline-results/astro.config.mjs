import baseConfig from "../../astro.config.base.mjs";
import icon from "astro-icon";
import { defineConfig } from "astro/config";
import pipelines_json from "./public/pipelines.json";

// Build redirects from pipelines data
let latestPipelineReleases = {};
pipelines_json.remote_workflows.map(
  (pipeline) => (latestPipelineReleases[pipeline.name] = `/${pipeline.name}/${pipeline.releases[0].tag_name}/`),
);

// https://astro.build/config
export default defineConfig({
  ...baseConfig,
  redirects: {
    ...latestPipelineReleases,
  },
  build: {
    ...baseConfig.build,
    assetsPrefix: "https://nf-core-pipeline-results.netlify.app/",
  },
  integrations: [
    ...baseConfig.integrations.filter((i) => i.name !== "astro-icon"),
    icon({
      iconDir: "../main-site/src/icons",
      include: baseConfig.integrations.find((i) => i.name === "astro-icon")?.include,
    }),
  ],
});
