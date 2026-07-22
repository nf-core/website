import type { APIRoute } from "astro";
import { octokit } from "@components/octokit";

export const GET: APIRoute = async () => {
    const versions = await octokit.paginate(octokit.rest.repos.listReleases, {
        owner: "nextflow-io",
        repo: "nextflow",
        per_page: 100,
    });

    const formattedVersions = versions
        .map((version) => ({
            version: version.tag_name,
            isEdge: version.prerelease,
            downloadUrl: version.assets[0]?.browser_download_url,
            downloadUrlAll: version.assets[1]?.browser_download_url,
            published_at: version.published_at,
        }))
        .sort((a, b) => {
            // remove `-edge` suffixes and then sort by nextflow's date-based versioning
            const aParts = a.version.replace("-edge", "").replace("v", "").split(".").map(Number);
            const bParts = b.version.replace("-edge", "").replace("v", "").split(".").map(Number);
            for (let i = 0; i < 3; i++) {
                if (aParts[i] !== bParts[i]) {
                    return bParts[i] - aParts[i];
                }
            }
            return 0;
        });

    // Fail if no data. setup-nextflow action depends on this file.
    if (formattedVersions.length === 0) {
        throw new Error("nextflow_version: GitHub returned no Nextflow versions; refusing to build an empty payload.");
    }

    return new Response(
        JSON.stringify({
            latest: {
                stable: formattedVersions.find((version) => !version.isEdge),
                edge: formattedVersions.find((version) => version.isEdge),
                everything: formattedVersions[0],
            },
            versions: formattedVersions,
        }),
        {
            status: 200,
            headers: {
                "Content-Type": "application/json",
            },
        },
    );
};
