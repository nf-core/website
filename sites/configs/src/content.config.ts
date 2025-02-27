import { z, defineCollection  } from 'astro:content';
import { githubFileLoader, type RenderedContent, md } from "@utils/loaders";

import type { AstroConfig } from "astro";


const configProcessor = async (text: string, config: AstroConfig): Promise<RenderedContent> => {
    let NFConfig = {
        config_profile_description: "",
        config_profile_contact: "",
        config_profile_url: "",
        executor: "",
    };

    const descMatch = text.match(/config_profile_description\s*=\s*("+|')([\s\S]*?)\1/);
    NFConfig.config_profile_description = descMatch?.[2] ?? "";

    const contactMatch = text.match(/config_profile_contact\s*=\s*(.*)/);
    NFConfig.config_profile_contact = contactMatch?.[1] ?? "";

    const urlMatch = text.match(/config_profile_url\s*=\s*(.*)/);
    NFConfig.config_profile_url = urlMatch?.[1] ?? "";

    const execMatch = text.match(/executor\s*=\s*(.*)/);
    NFConfig.executor = execMatch?.[1] ?? "";
    if (NFConfig.executor === "") {
        const execNameMatch = text.match(/executor\s*{\s*name\s*=\s*(.*)\s*/);
        NFConfig.executor = execNameMatch?.[1] ?? "";
    }
    if (NFConfig.executor?.includes("?") && NFConfig.executor?.includes(":")) {
        const match = NFConfig.executor.match(/{.*\?(.*):(.*)}/);
        NFConfig.executor = match
            ? match
                  .slice(1, 3)
                  .map((item) => item.trim().replace(/'/g, ""))
                  .join(", ")
            : "";
        // replace trailing comments (starting with '//')
        NFConfig.executor = NFConfig.executor.replace(/\/\/.*$/, "");
    }
    //remove double and quotes
    Object.keys(NFConfig).forEach((key) => {
        NFConfig[key] = NFConfig[key]?.replace(/'|"/g, "");
    });

    //remove provided by nf-core/configs string from description
    NFConfig.config_profile_description = NFConfig.config_profile_description.replace(" provided by nf-core/configs", "");

    return {
        html: text,
        metadata: NFConfig,
    };
};



const configs = defineCollection({
    loader: githubFileLoader({
        org: "nf-core",
        repo: "configs",
        ref: "master",
        path: /\.(md|mdx|config)$/,
        processors: {
            md,
            config: configProcessor,
        },
        getDates: true,
    }),
    schema: z.object({
        extension: z.literal('md').or(z.literal('config')),
        sha: z.string(),
        lastCommit: z.string().datetime(),
    }),
});

export const collections = {
    configs,
};
