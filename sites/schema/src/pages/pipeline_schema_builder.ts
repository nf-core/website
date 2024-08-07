import type { APIRoute } from 'astro';
export const prerender = false;

export const POST: APIRoute = async ({ request }) => {
    console.log('Received request');

    console.log('Received JSON request');
    console.log('Request:', request);
    const body = await request.json();
    console.log('Body:', body);
    const schema = body.schema.schema;
    return new Response(
        JSON.stringify({
            message: 'Your schema was: ' + schema,
        }),
        {
            status: 200,
        },
    );
    return new Response(null, { status: 400 });
};
