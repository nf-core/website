# AGENTS.md

## Cursor Cloud specific instructions

This repo is the **nf-core website**: an [Astro](https://astro.build/) static-site
mono-repo split into independent sub-sites under `sites/*` wired together with npm
workspaces. The update script runs `npm install` on startup (installs both the root
and all workspace dependencies), so the notes below focus on running/testing.

### Services (Astro sub-sites)

Each `sites/<name>` is its own Astro app with the same scripts (`dev`, `build`,
`test`). Run one at a time; they are disjoint, so cross-site links (e.g. opening
`/pipelines` while running `main-site`) will 404 in dev — this is expected.

- `sites/main-site` – main website (blog, events, components). Works with no token.
- `sites/docs` – documentation. Works with no token.
- `sites/pipelines`, `sites/modules-subworkflows`, `sites/configs`,
  `sites/pipeline-results`, `sites/stickers` – pull data from the GitHub API.
  Pre-built JSON in `public/` lets them run for basic dev, but to rebuild data or
  avoid API rate limits add a `GITHUB_TOKEN` to a root `.env` and symlink it into
  the sub-site, e.g. `ln -s .env sites/pipelines/.env` (see `README.md`).

### Run (dev)

`npm run dev --workspace sites/main-site` serves at http://localhost:4321/.
Add `-- --host` to expose it on the network. First start is slow (~10s) because
Astro syncs content and optimizes deps.

### Lint

The CI lint gate is **prettier driven by prek** (the `@j178/prek` binary), only over
`*.astro,*.svelte,*.mdx,*.md,*.yml,*.yaml`. Run it with
`./node_modules/.bin/prek run --all-files`.
Note: the per-workspace `npm run lint` (`prettier --check .`) scans _all_ files
including vendored FontAwesome CSS and reports long-standing warnings that are NOT
part of the CI gate — prefer `prek` to judge lint status.

### Test

Playwright lives in `sites/main-site`: `npm run test -w sites/main-site`
(add `-- --project=chromium` to run a single browser). Browsers must be installed
once with `npx playwright install --with-deps chromium`. Locally the tests hit the
dev server at `http://localhost:4321`, so start `main-site` dev first. CI instead
points `PLAYWRIGHT_TEST_BASE_URL` at a Netlify deploy preview.

### Gotcha

Root devDependencies (`prek`, `prettier`) are only installed by plain `npm install`,
**not** `npm install --workspaces` (which the README mentions). The update script
uses plain `npm install` for this reason. `npm install` may add harmless `"peer": true`
annotations to `package-lock.json`; leave those uncommitted.
