#! /usr/bin/env node
// move markdown files from cache to src/content

import { promises, readdirSync, statSync, existsSync, mkdirSync } from 'fs';

import path from 'path';
export async function buildContentCollection(){
// go through all files in .cache, move them to src/content

  const getAllMDFiles = dir =>
    readdirSync(dir).reduce((files, file) => {
        const name = path.join(dir, file);
        const isDirectory = statSync(name).isDirectory();
        if (isDirectory){
            return [...files, ...getAllMDFiles(name)];
        } else {
          if (/\.mdx?$/.test(name)){
            return [...files, name];
          }
          return files;
        }

        // check if file is markdown

        // return isDirectory ? [...files, ...getAllMDFiles(name)] : [...files, /\.mdx?$/.test(name)? name: ];
    }, []);

  const files = getAllMDFiles('.cache');
  if (!existsSync('src/content/pipelines')){
      mkdirSync('src/content/pipelines', { recursive: true });
    }
  Promise.all(
    // create src/content/pipelines folder if it doesn't exist


    files.map(async (f) => {
      let content = await promises.readFile(f, 'utf8');
      const pathParts = f.split('/');
      const pipeline = pathParts[1];
      const version = pathParts[2];
      // make relative links to other markdown files absolute to current Astro.url.pathname
        content = content.replaceAll(/\[([^\]]+)\]\((?!http)(?!#)(.*?)\)/g, (match, p1, p2) => {
            const link = path.join(pipeline, version, p2);
            return `[${p1}](/${link})`;
        });

      const newPath = f.replace('.cache', 'src/content/pipelines');
      const parent = newPath.split('/').slice(0, -1).join('/');
      await promises.mkdir(parent, { recursive: true });
      await promises.writeFile(newPath, content);
    }),
  );

};

buildContentCollection();
