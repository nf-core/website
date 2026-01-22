---
title: Manually syncing your pipeline
subtitle: How to manually synchronise your nf-core pipeline with the latest nf-core/tools template
shortTitle: Manually syncing your pipeline
weight: 3
---

In rare cases, you need to trigger synchronisation manually. This applies when automated sync did not run during an `nf-core/tools` release on GitHub, when you want to perform a targeted sync, or when you maintain a custom pipeline created with `nf-core pipelines create`.

:::note
Automated template syncs only apply to official nf-core pipelines. Custom pipelines require manual synchronisation.
:::

## Sync an official nf-core pipeline

Use this approach when working with an official nf-core pipeline that needs manual synchronisation:

1. Navigate to your pipeline directory and check out your development branch:

   ```bash
   cd <pipeline_name>
   git checkout dev
   ```

1. If you work on a fork, pull the `TEMPLATE` branch from the upstream nf-core repository:

   ```bash
   git remote add upstream https://github.com/nf-core/<pipeline_name>.git
   git checkout --track upstream/TEMPLATE
   git checkout dev
   ```

1. Run the sync command:

   ```bash
   nf-core pipelines sync
   ```

The sync command creates a merge that you can resolve locally, similar to automated sync pull requests.

## Sync a custom pipeline with a pull request

Use this approach when you maintain a custom pipeline and want to create a pull request to your `dev` branch:

1. Navigate to your pipeline directory:

   ```bash
   cd <pipeline_name>
   ```

1. Run the sync command with pull request options:

   ```bash
   nf-core pipelines sync \
      --dir <pipeline_directory> \
      --from-branch dev \
      --pull-request \
      --username <github_username> \
      --github-repository <github_repository_url>
   ```

This creates a pull request in your repository that you can review and merge through the GitHub interface.
