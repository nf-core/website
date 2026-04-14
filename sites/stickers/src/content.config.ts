import { z, defineCollection } from "astro:content";
import type { AstroConfig } from "astro";
import { githubFileLoader } from "@utils/loaders";
import { octokit } from "@components/octokit.js";

// Fetch teams data for team sticker descriptions
const teams = await octokit.request("GET /orgs/{org}/teams", {
  org: "nf-core",
  headers: {
    "X-GitHub-Api-Version": "2022-11-28",
  },
});

const stickers = defineCollection({
  loader: githubFileLoader({
    org: "nf-core",
    repo: "logos",
    ref: "master",
    path: /^hexagon-stickers\/.+\.png$/,
    getDates: false,
    processors: {
      png: async (text: string, config: AstroConfig, filepath?: string) => {
        if (!filepath) {
          console.log("filepath is undefined");
          return {
            html: "",
            metadata: {},
          };
        }
        // skip template and preview images
        if (filepath.includes("-TEMPLATE-") || filepath.includes("-preview")) {
          console.log("skipping template or preview image", filepath);
          return {
            html: "",
            metadata: {},
          };
        }
        // Extract category and name from the file path
        const pathMatch = filepath.match(/^hexagon-stickers\/([^/]+)(?:\/([^/]+))?\/[^/]+\.png$/);
        if (!pathMatch) {
          return {
            html: "",
            metadata: {
              name: "unknown",
              category: "other",
              imageUrl: "",
              alt: "Unknown sticker",
              link: "https://nf-co.re",
              description: "Unknown sticker",
            },
          };
        }

        const [, category, subdir] = pathMatch;
        const name = subdir || category;
        const imageUrl = `https://raw.githubusercontent.com/nf-core/logos/master/${filepath}`;
        let description = subdir ? `${name} pipeline sticker` : "";
        if (category === "pipelines") {
          const response = await octokit.request("GET /repos/{owner}/{repo}", {
            owner: "nf-core",
            repo: name,
            headers: {
              "X-GitHub-Api-Version": "2022-11-28",
            },
          });
          description = response.data.description || description;
        }
        if (category === "teams") {
          description =
            teams.data.find(
              (team) => team.name === filepath.split("/").pop()?.replace("nf-core-", "").replace(".png", ""),
            )?.description || description;
        }

        return {
          html: "",
          metadata: {
            name,
            category,
            imageUrl,
            alt: `${name} sticker`,
            link: subdir ? `https://github.com/nf-core/${name}` : "https://nf-co.re",
            description,
          },
        };
      },
    },
  }),
  schema: z.object({
    name: z.string(),
    category: z.string(),
    imageUrl: z.string(),
    alt: z.string(),
    link: z.string().url(),
    description: z.string(),
  }),
});

export const collections = {
  stickers,
};
