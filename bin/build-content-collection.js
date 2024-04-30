#! /usr/bin/env node
// move markdown files from cache to src/content
import { promises, readdirSync, statSync, existsSync, mkdirSync } from 'fs';
import path from 'path';


(async () => {
async function buildContentCollection() {
  // go through all files in .cache, move them to src/content
  console.log('Building content collection');

  const getAllMDFiles = (dir) =>
    readdirSync(dir).reduce((files, file) => {
      const name = path.join(dir, file);
      const isDirectory = statSync(name).isDirectory();
      if (isDirectory) {
        return [...files, ...getAllMDFiles(name)];
      } else {
        if (/\.mdx?$/.test(name)) {
          return [...files, name];
        }
        return files;
      }
    }, []);

      console.log('current path:', path.resolve());
      const files = getAllMDFiles('.cache');
      if (!existsSync('./sites/pipelines/src/content/pipelines')) {
        mkdirSync('./sites/pipelines/src/content/pipelines', { recursive: true });
        console.log('Created sites/pipelines/src/content/pipelines folder');
        console.log('current path:', path.resolve());
      }
  Promise.all(
    // create sites/pipelines/src/content/pipelines folder if it doesn't exist

    files.map(async (f) => {
      let content = await promises.readFile(f, 'utf8');
      const pathParts = f.split('/');
      const pipeline = pathParts[1];
      const version = pathParts[2];
      // make relative links to png and svg files absolute in markdown to current Astro.url.pathname
      content = content.replaceAll(
        /(\]\()(docs\/images\/.*?\.png|svg)/gmu,
        `$1${`https://raw.githubusercontent.com/nf-core/${pipeline}/${version}/$2`}`
      );
      const newPath = f.replace('.cache', 'sites/pipelines/src/content/pipelines');
      const parent = newPath.split('/').slice(0, -1).join('/');
      await promises.mkdir(parent, { recursive: true });
      await promises.writeFile(newPath, content);
    })
  );
  return true;
};

await buildContentCollection();
})();
