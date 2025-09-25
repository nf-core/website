import test from 'ava';
import { execa } from 'execa';
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const mdxFolder = path.join(__dirname, './mdx');

test('Convert MDX to HTML', async t => {
    const cases = fs.readdirSync(mdxFolder);

    for (const caseFolder of cases) {
        const casePath = path.join(mdxFolder, caseFolder);
        if (fs.statSync(casePath).isDirectory()) {
            const mdxFilePath = path.join(casePath, 'test.mdx');
            const htmlFilePath = path.join(casePath, 'test.html');

            const expectedHtmlContent = fs.readFileSync(htmlFilePath, 'utf8').replace(/\s+/g, '');;
            const { stdout: actualHtmlContent } = await execa('./bin/cli.js', [mdxFilePath]);

            t.is(actualHtmlContent.replace(/\s+/g, ''), expectedHtmlContent, `Failed for case: ${caseFolder}`);
        }
    }
});
