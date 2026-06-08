#! /usr/bin/env node
// Checks the main site's machine-discovery catalog against what's actually in public/.
//
// The api-catalog (sites/main-site/public/.well-known/api-catalog, an RFC 9727 link set) is
// the source of truth for machine-readable resources. Three checks, all static, no build:
//   1. Drift    — every static file the catalog links to actually exists.
//   2. Coverage — every static discovery file in public/ is listed in the catalog, so a
//                 newly-added .json can't end up silently undiscoverable.
//   3. CORS     — every static discovery resource has an Access-Control-Allow-Origin rule in
//                 netlify.toml (else browser-based agents can't read it cross-origin), and no
//                 stale discovery CORS rules linger. Same resource list, kept in two files.
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
// The catalog's own path — a discovery resource, but not listed inside its own linkset.
const CATALOG_PATH = '/.well-known/api-catalog';
// netlify.toml CORS rules that aren't discovery resources: the asset modules, and the
// homepage (whose ACAO exists so cross-origin JS can read its Link header).
const CORS_INFRA = new Set(['/_astro/*', '/']);

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

// 3. CORS: the catalog's static resources and netlify.toml's CORS rules are the same list.
const needsCors = new Set([CATALOG_PATH, ...[...linked].filter(isStatic)]);
const corsPaths = new Set();
const tomlCode = readFileSync(NETLIFY, 'utf8')
  .split('\n')
  .filter((l) => !l.trim().startsWith('#')) // drop comments (they mention the header name)
  .join('\n');
for (const block of tomlCode.split('[[headers]]').slice(1)) {
  const forPath = block.match(/for\s*=\s*"([^"]+)"/)?.[1];
  if (forPath && /Access-Control-Allow-Origin/.test(block)) corsPaths.add(forPath);
}
for (const r of needsCors) {
  if (!corsPaths.has(r)) {
    errors.push(`${r} is a discovery resource but has no Access-Control-Allow-Origin rule in netlify.toml.`);
  }
}
for (const p of corsPaths) {
  if (!CORS_INFRA.has(p) && !needsCors.has(p)) {
    errors.push(
      `netlify.toml grants CORS to ${p}, which isn't a discovery resource in the api-catalog ` +
        `(remove the rule, add the resource to the catalog, or list it in CORS_INFRA).`,
    );
  }
}

if (errors.length) {
  console.error('✗ api-catalog check failed:');
  for (const e of errors) console.error('  • ' + e);
  process.exit(1);
}
console.log('✓ api-catalog matches the discovery files in public/ and netlify.toml CORS rules.');
