import baseConfig from "../../astro.config.base.mjs";
import icon from "astro-icon";
import { defineConfig } from "astro/config";
import pipelines_json from "./public/pipelines.json";

// Build redirects from pipelines data
let latestPipelineReleases = {};
pipelines_json.remote_workflows.map(
  (pipeline) => (latestPipelineReleases[pipeline.name] = `/${pipeline.name}/${pipeline.releases[0].tag_name}/`),
);

let pipelineResults = {};
pipelines_json.remote_workflows.map(
  (pipeline) =>
    (pipelineResults[`/${pipeline.name}/:version/results/*`] =
      `https://nf-core-pipeline-results.netlify.app/${pipeline.name}/:version/results/:splat 200!`),
);

// https://astro.build/config
export default defineConfig({
  ...baseConfig,
  redirects: {
    ...latestPipelineReleases,
    ...pipelineResults,
  },
  build: {
    ...baseConfig.build,
    assetsPrefix: "https://nf-core-pipelines.netlify.app/",
  },
  integrations: [
    ...baseConfig.integrations.filter((i) => i.name !== "astro-icon"),
    icon({
      iconDir: "../main-site/src/icons",
      include: baseConfig.integrations.find((i) => i.name === "astro-icon")?.include,
    }),
  ],
});
