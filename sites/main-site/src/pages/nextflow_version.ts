import type { APIRoute } from 'astro';
import octokit from '@components/octokit';
import type { OctokitResponse } from '@octokit/types';

export const GET: APIRoute = async ({ params, request }) => {
    try {
        // get all releases as versions, go through all release pages

        const versions: {}[] = [];
        let page = 1;
        let releases: OctokitResponse<any>;
        do {
            releases = await octokit.request('GET /repos/{owner}/{repo}/releases', {
                owner: 'nextflow-io',
                repo: 'nextflow',
                page,
            });
            versions.push(...releases.data);
            page++;
        } while (releases.headers.link?.includes('rel="next"'));

        const formattedVersions = versions.map((version: any) => ({
            version: version['tag_name'],
            isEdge: version['prerelease'],
            downloadUrl: version['assets'][0] && version['assets'][0]['browser_download_url'],
            downloadUrlAll: version['assets'][1] && version['assets'][1]['browser_download_url'],
            published_at: version['published_at'],
        })).sort((a: any, b: any) => { // sort version by version number, by removing "-edge" suffices and then comparing semver
            const aVersion = a['version'].replace('-edge', '').replace('v', '');
            const bVersion = b['version'].replace('-edge', '').replace('v', '');
            // check first major, then minor, then patch
            const aVersionSplit = aVersion.split('.');
            const bVersionSplit = bVersion.split('.');
            for (let i = 0; i < 3; i++) {
                if (parseInt(aVersionSplit[i]) > parseInt(bVersionSplit[i])) {
                    return -1;
                } else if (parseInt(aVersionSplit[i]) < parseInt(bVersionSplit[i])) {
                    return 1;
                }
            }
            return 0;
        });

        return new Response(
            JSON.stringify({
                latest:{
                    stable: formattedVersions.find((version: any) => !version.isEdge),
                    edge: formattedVersions.find((version: any) => version.isEdge),
                    everything: formattedVersions[0],
                },
                versions: formattedVersions,
            }), {
        status: 200,
        headers: {
        "Content-Type": "application/json"
        }
    }
        );
    } catch (error) {
        return new Response(
            JSON.stringify({
                error: error.message,
            }),
            { status: 500 },
        );
    }
};
