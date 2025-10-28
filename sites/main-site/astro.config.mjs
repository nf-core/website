import baseConfig from "../../astro.config.base.mjs";
import icon from "astro-icon";
import { defineConfig } from "astro/config";
import pipelines_json from "./public/pipelines.json";

// Build redirects from pipelines data
let pipelineRedirects = {};
pipelines_json.remote_workflows.map((pipeline) => {
  pipelineRedirects[`/${pipeline.name}/:version/*`] =
    `https://nf-core-pipelines.netlify.app/${pipeline.name}/:version/:splat 200!`;
  pipelineRedirects[`/${pipeline.name}/`] = `https://nf-core-pipelines.netlify.app/${pipeline.name} 200!`;
});

// https://astro.build/config
export default defineConfig({
  ...baseConfig,
  redirects: {
    ...pipelineRedirects,
  },
  build: {
    ...baseConfig.build,
    assetsPrefix:
      process.env.CONTEXT === "production" ? "https://nf-core-main-site.netlify.app/" : process.env.DEPLOY_PRIME_URL,
  },
  vite: {
    ...baseConfig.vite,
    envPrefix: ["PUBLIC_"], // Remove GITHUB_ prefix to prevent token exposure
  },
  image: {
    ...baseConfig.image,
    domains: [...baseConfig.image.domains, "netlify.app", "nf-co.re", "nf-core-main-site.netlify.app"],
  },
  integrations: [
    ...baseConfig.integrations.filter((i) => i.name !== "astro-icon"),
    icon({
      include: {
        "file-icons": ["nextflow"],
        logos: [
          "bluesky",
          "twitter",
          "mastodon-icon",
          "slack-icon",
          "aws",
          "microsoft-azure",
          "github-actions",
          "youtube-icon",
          "linkedin",
        ],
        fa: ["github"],
        "fa-brands": ["github"],
        mdi: [
          "aws",
          "linkedin",
          "slack",
          "youtube",
          "cloud-outline",
          "timeline-check-outline",
          "book-information-variant",
          "package-variant",
          "progress-check",
        ],
        octicon: [
          "chevron-right-16",
          "git-pull-request-16",
          "law-16",
          "link-external-16",
          "mortar-board-16",
          "play-16",
          "table-16",
          "tasklist-16",
          "terminal-16",
          "tools-16",
        ],
        "simple-icons": ["bluesky"],
        ri: ["open-source-line"],
        "pepicons-print": ["t-shirt"],
        fluent: ["paint-brush-sparkle-24-filled"],
        bi: ["bag-heart-fill"],
      },
    }),
  ],
});
