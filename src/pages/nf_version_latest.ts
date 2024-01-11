import type { APIRoute } from 'astro';


export const GET: APIRoute = async ({ params, request }) => {
    try {
        // get latest stable nextflow version from github releases
        const response = await fetch('https://api.github.com/repos/nextflow-io/nextflow/releases/latest');
        if (!response.ok) {
            throw new Error('Failed to fetch latest version');
        }
        const data = await response.json();
        return new Response(
            JSON.stringify({
                'version': data['tag_name']
            })
        );
    } catch (error) {
        return new Response(
            JSON.stringify({
                'error': error.message
            }),
            { status: 500 }
        );
    }
}
