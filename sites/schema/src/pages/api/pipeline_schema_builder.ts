import type { APIRoute } from 'astro';

interface SchemaData {
    post_content: string;
    api: string;
    version: string;
    status: string;
    schema: string;
}

interface ProcessedSchemaData extends Omit<SchemaData, 'schema'> {
    schema: object;
}

export const post: APIRoute = async ({ request }) => {
    const data: SchemaData = await request.json();

    const processedData: ProcessedSchemaData = {
        post_content: data.post_content,
        api: data.api,
        version: data.version,
        status: data.status,
        schema: JSON.parse(data.schema),
    };

    return new Response(JSON.stringify(processedData), {
        status: 200,
        headers: {
            'Content-Type': 'application/json',
        },
    });
};
