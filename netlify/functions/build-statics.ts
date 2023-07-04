import { Handler } from 'aws-lambda';
import { exec } from 'child_process';

export const handler: Handler = async () => {
    try {
        console.log('Building cache...');
        await runCommand('npm run build-pipeline-json');
        await runCommand('npm run build-component-json');
        await runCommand('node bin/build-cache.js');
        await runCommand('tar -cJf --no-xattrs -f .cache.tar.xz .cache');

        return {
            statusCode: 200,
            body: 'Commands executed successfully.',
        };
    } catch (error) {
        return {
            statusCode: 500,
            body: error.message,
        };
    }
};

function runCommand(command: string): Promise<void> {
    return new Promise<void>((resolve, reject) => {
        exec(command, (error, stdout, stderr) => {
            if (error) {
                reject(new Error(`Command "${command}" failed: ${error.message}`));
            } else {
                resolve();
            }
        });
    });
}
