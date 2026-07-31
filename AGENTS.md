# AGENTS.md

This file provides guidance to AI coding agents (Claude Code and others) when working with code in this repository.

The source for **https://nf-co.re** — the nf-core community website. Built with **Astro** (static site generation), **Svelte** (interactive components) and **Bootstrap** (CSS), managed as **npm workspaces** with **Node 22**.

## Architecture

The repo is a **monorepo of independent Astro sites**, one per subdirectory of `sites/*`, each an npm workspace ([background blog post](https://nf-co.re/blog/2024/new-website-structure)):

- `main-site` — the apex `nf-co.re`: homepage, blog, events, community/governance pages, and shared components/layouts that other sites import.
- `docs` — documentation.
- `pipelines`, `modules-subworkflows`, `pipeline-results`, `configs`, `stickers` — data-driven sites that pull content from the **GitHub API** at build time.

You usually work on **one sub-site at a time**, not the whole repo.

Two cross-cutting things are non-obvious and worth knowing before editing build config:

- **Each sub-site deploys as its own separate Netlify site**, and the apex `main-site` stitches them together: it **proxies** many URL paths to the other sites via rewrite rules in `sites/main-site/public/_redirects` (e.g. `/usage/*` and `/contributing/*` are served by `docs`).
- **Assets load cross-origin.** Each sub-site sets `build.assetsPrefix` (in its `astro.config.mjs`) to its own Netlify subdomain, so a page served under `nf-co.re` loads its JS/CSS/fonts from a different origin. This is why each sub-site's `netlify.toml` serves `/_astro/*` with an `Access-Control-Allow-Origin` CORS header.

Shared build logic lives in `bin/` (Astro/remark/rehype plugins, and the scripts that generate the machine-readable `public/pipelines.json` and `public/components.json` catalogues). Common Astro config is in `astro.config.base.mjs`, extended by each sub-site.

## Development

```sh
npm install --workspaces                 # install dependencies for all sub-sites
npm run dev --workspace sites/main-site  # dev server for one sub-site (astro dev)
npm run build --workspace sites/docs     # production build for one sub-site
npm run test-all                         # Playwright tests across all sub-sites
```

Tests are **Playwright** (`npm run test --workspace sites/<name>`, or `npx playwright test <file>` within a sub-site for a single test).

**GitHub token:** the data-driven sub-sites (`pipelines`, `pipeline-results`, `configs`, `modules-subworkflows`, `stickers`) hit the GitHub API during the build. Add a `.env` at the repo root with a `GITHUB_TOKEN` (a personal access token with `public_repo` scope) to avoid rate limits.

## Content & contribution guides

See [`README.md`](README.md) for step-by-step instructions and examples on adding content, including: adding a blog post, event, security advisory, or announcement banner; updating the JSON files; adding a new sub-site to the monorepo; the file structure; and the blog post style/contribution guidelines.
