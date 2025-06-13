import type { APIRoute } from 'astro';
import { getNextflowVersions } from '@utils/functions';

export const GET: APIRoute = async () => {
    try {

        const nextflowVersions = await getNextflowVersions();

        return new Response(
            JSON.stringify({
                latest: {
                    stable: nextflowVersions.find(version => !version.isEdge),
                    edge: nextflowVersions.find(version => version.isEdge),
                    everything: nextflowVersions[0],
                },
                versions: nextflowVersions,
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
