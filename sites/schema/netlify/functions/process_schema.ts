import type { Handler, HandlerEvent, HandlerContext } from '@netlify/functions';
import { getStore, type Store } from '@netlify/blobs';

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
    if (event.httpMethod === 'POST') {
        if (!event.body) {
            return { statusCode: 400, body: 'Bad Request' };
        }
        if (event.headers['content-type'] !== 'application/json') {
            const formData = new URLSearchParams(event.body);
            const data: SchemaData = {
                post_content: formData.get('post_content') || '',
                api: formData.get('api') || '',
                version: formData.get('version') || '',
                status: formData.get('status') || '',
                schema: JSON.parse(formData.get('schema') || ''),
            };

            // Store the processed data in a Netlify Blob
            const store = getStore({
                name: 'schema',
                siteID: process.env.SITE_ID,
                token: process.env.NETLIFY_AUTH_TOKEN,
            });
            const key = `schema_${Date.now()}_${Math.random().toString().substring(3)}`; // Use a unique key, based on timestamp and a random string
            await store.setJSON(key, data);

            return {
                statusCode: 200,
                body: JSON.stringify({
                    message: 'Data stored successfully',
                    key,
                    status: 'received',
                    web_url: 'http://localhost:8888/schema_builder?id=' + key,
                    api_url: 'http://localhost:8888/.netlify/functions/process_schema?id=' + key,
                }),
            };
        } else {
            const data: SchemaData = JSON.parse(event.body);
            const key = event.queryStringParameters?.id;
            // Store the processed data in a Netlify Blob
            const store = getStore({
                name: 'schema',
                siteID: process.env.SITE_ID,
                token: process.env.NETLIFY_AUTH_TOKEN,
            });
            if (key) {
                await store.setJSON(key, data);
            }
        }
    }
    if (event.httpMethod === 'GET') {
        const key = event.queryStringParameters?.id;
        if (!key) {
            return { statusCode: 400, body: 'Bad Request' };
        }
        const store = getStore({ name: 'schema', siteID: process.env.SITE_ID, token: process.env.NETLIFY_AUTH_TOKEN });

        // Retrieve the processed data from the Netlify Blob
        const schema = await store.get(key);
        if (!schema) {
            return { statusCode: 404, body: 'Not Found' };
        }
        const data = {
            schema: JSON.parse(schema),
        };

        if (data.schema.status === 'waiting_for_user') {
            return {
                statusCode: 200,
                body: JSON.stringify({
                    message: 'GET request received',
                    status: 'waiting_for_user',
                    data,
                }),
            };
        } else if (data.schema.status === 'processed') {
            return {
                statusCode: 200,
                body: JSON.stringify({
                    message: 'GET request received',
                    status: 'processed',
                    data,
                }),
            };
        }
        return { statusCode: 400, body: 'Bad Request' };
    }
};

export { handler };
