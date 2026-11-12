#! /usr/bin/env node
// Checks the main site's machine-discovery catalog against what's actually in public/.
//
// The api-catalog (sites/main-site/public/.well-known/api-catalog, an RFC 9727 link set) is
// the source of truth for machine-readable resources. Three checks, all static, no build:
//   1. Drift    — every static file the catalog links to actually exists.
//   2. Coverage — every static discovery file in public/ is listed in the catalog, so a
//                 newly-added .json can't end up silently undiscoverable.
//   3. CORS     — the api-catalog has an Access-Control-Allow-Origin rule in netlify.toml, so
//                 a browser-based agent can fetch it cross-origin. The catalog is the single
//                 entry point; the resources it lists are served same-origin everywhere and
//                 read server-side by crawlers, so they don't need their own CORS rules.
//
// Only .json/.txt files are checked, since those are the discovery resources that exist as
// static source. Page routes (/docs/) and build-generated files (sitemap-index.xml) aren't
// static source, so they're out of scope — verifying those resolve is the live
// broken-link-checker's job against the Netlify preview.

import { readFileSync, existsSync, readdirSync } from 'fs';
import path from 'path';

const SITE = 'https://nf-co.re';
const PUBLIC = path.join(path.resolve(), 'sites/main-site/public');
const CATALOG = path.join(PUBLIC, '.well-known/api-catalog');
const NETLIFY = path.join(PUBLIC, '..', 'netlify.toml');

// Extensions that exist as static source files in public/, so can be checked on disk.
const STATIC_EXTS = new Set(['.json', '.txt']);
// public/ files that needn't appear in the catalog (not machine-readable resources).
const IGNORE = new Set(['robots.txt']);
// The catalog's served path — must be cross-origin readable so agents can fetch it.
const CATALOG_PATH = '/.well-known/api-catalog';

const isStatic = (p) => STATIC_EXTS.has(path.extname(p));
const errors = [];

if (!existsSync(CATALOG)) {
  console.error(`✗ api-catalog not found at ${path.relative('.', CATALOG)}`);
  process.exit(1);
}

// Collect every on-site path the catalog links to, as site-relative "/path" strings.
const linked = new Set();
const catalog = JSON.parse(readFileSync(CATALOG, 'utf8'));
for (const entry of catalog.linkset ?? []) {
  for (const [rel, links] of Object.entries(entry)) {
    if (rel === 'anchor') continue; // the anchor is the base URL, not a link
    for (const { href } of links) {
      if (href?.startsWith(SITE + '/')) linked.add(href.slice(SITE.length));
    }
  }
}

// 1. Drift: catalog links a static file that doesn't exist.
for (const link of linked) {
  if (isStatic(link) && !existsSync(path.join(PUBLIC, link))) {
    errors.push(`Catalog links ${SITE}${link} but public${link} does not exist.`);
  }
}

// 2. Coverage: a static file exists in public/ but the catalog doesn't list it.
for (const name of readdirSync(PUBLIC)) {
  if (isStatic(name) && !IGNORE.has(name) && !linked.has('/' + name)) {
    errors.push(
      `public/${name} looks machine-readable but isn't in the api-catalog ` +
        `(add it to the catalog, or to IGNORE in bin/check-api-catalog.js).`,
    );
  }
}

// 3. CORS: the api-catalog must be cross-origin readable so browser-based agents can fetch it.
const corsPaths = new Set();
const tomlCode = readFileSync(NETLIFY, 'utf8')
  .split('\n')
  .filter((l) => !l.trim().startsWith('#')) // drop comments (they mention the header name)
  .join('\n');
for (const block of tomlCode.split('[[headers]]').slice(1)) {
  const forPath = block.match(/for\s*=\s*"([^"]+)"/)?.[1];
  if (forPath && /Access-Control-Allow-Origin/.test(block)) corsPaths.add(forPath);
}
if (!corsPaths.has(CATALOG_PATH)) {
  errors.push(`${CATALOG_PATH} has no Access-Control-Allow-Origin rule in netlify.toml; agents can't fetch it cross-origin.`);
}

if (errors.length) {
  console.error('✗ api-catalog check failed:');
  for (const e of errors) console.error('  • ' + e);
  process.exit(1);
}
console.log('✓ api-catalog matches the discovery files in public/ and netlify.toml CORS rules.');
