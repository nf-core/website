import type { Handler, HandlerEvent, HandlerContext } from '@netlify/functions';
import { getStore } from '@netlify/blobs';

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

const handler: Handler = async (event: HandlerEvent, context: HandlerContext) => {
    if (event.httpMethod !== 'POST') {
        return { statusCode: 405, body: 'Method Not Allowed' };
    }
    if (!event.body) {
        return { statusCode: 400, body: 'Bad Request' };
    }
    const formData = new URLSearchParams(event.body);
    const data: SchemaData = {
        post_content: formData.get('post_content') || '',
        api: formData.get('api') || '',
        version: formData.get('version') || '',
        status: formData.get('status') || '',
        schema: formData.get('schema') || '',
    };

    const processedData: ProcessedSchemaData = {
        post_content: data.post_content,
        api: data.api,
        version: data.version,
        status: data.status,
        schema: JSON.parse(data.schema),
    };

    // Store the processed data in a Netlify Blob
    const store = getStore('schema');
    const key = `schema_${Date.now()}`; // Use a unique key, e.g., based on timestamp
    await store.set(key, JSON.stringify(processedData));

    return {
        statusCode: 200,
        body: JSON.stringify({ message: 'Data stored successfully', key, status: 'received' }),
    };
};

export { handler };
