---
title: Synchronisation
subtitle: Learn how nf-core pipelines are automatically kept up to date with community standards.
---

# Introduction

To keep all the nf-core pipelines up-to-date with the latest version of the community standards, we have implemented a synchronisation tool.
This ensures that updates to the community standards are propagated to all nf-core pipelines.

There are three topics covered in this documentation page:

1. [Merging automated PRs](#merging-automated-prs)
2. [Manual synchronisation](#manual-synchronisation)
3. [Setting up a pipeline for syncing retrospectively](#setting-up-a-pipeline-for-syncing-retrospectively)
4. [Fixing a broken `TEMPLATE` branch](#fixing-a-broken-template-branch)

### How template synchronisation works

The `nf-core` helper tools have a subcommand for synchronising a pipeline with the nf-core template (`nf-core sync`).
Although this can be run manually, it is usually only used by the GitHub Actions automation:
when a new version of `nf-core/tools` is released it runs for all nf-core pipelines and automatically opens pull-requests (PRs) with the necessary changes required to update every pipeline.
These pull requests then need to be manually resolved and merged by the pipeline maintainers.

Behind the scenes, this synchronisation is done by using `git`.
Each repository has a special `TEMPLATE` branch which contains only the "vanilla" code made by the `nf-core create` tool.
The synchronisation tool fetches the essential variables needed to recreate the pipeline and uses this to trigger a `nf-core create --no-git` command with the latest version of the template.
The result from this is then compared against what is stored in the `TEMPLATE` branch and committed. During an automated sync, a copy of the `TEMPLATE` branch called `nf-core-template-merge-<version>` will be made (to avoid `dev` history to end up in `TEMPLATE` branch after solving merge conflicts), and a PR from this new branch will be opened against your `dev`.
When merging from the `nf-core-template-merge-<version>` branch back into the main `dev` branch of the pipeline, `git` should be clever enough to know what has changed since the template was first used, and therefore, it will only present the relevant changes.

For this to work in practice, the `TEMPLATE` branch needs to have a shared `git` history with the `master` branch of the pipeline.
The `nf-core create` command initially does this by enforcing a first commit to the `master` branch before any development has taken place.
If the pipeline _was not_ created by the `nf-core create` command, this has to be set up manually.
For instructions on this, see [Setting up a pipeline for syncing retrospectively](#setting-up-a-pipeline-for-syncing-retrospectively).
