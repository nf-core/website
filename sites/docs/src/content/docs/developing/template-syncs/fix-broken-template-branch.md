---
title: Fixing a broken TEMPLATE branch
description: How to fix a broken TEMPLATE branch
shortTitle: Fixing a broken TEMPLATE branch
---

If you resolve merge conflicts through the GitHub web interface, the commit history from your `dev` branch merges into `TEMPLATE`.
This causes problems in future template syncs because the updated `TEMPLATE` branch removes pipeline-specific files that were accidentally included.
You may see over 100 files with conflicts during the next sync.

## Identify a broken TEMPLATE branch

Check whether `dev` was accidentally merged into `TEMPLATE`:

1. Review the number of commits in your `TEMPLATE` branch.

:::note
A healthy `TEMPLATE` branch typically has approximately 5-15 commits, depending on how many template updates your pipeline has received.
:::

1. Check your repository's **Network Graph** under the **Insights** tab.
1. Look through the `TEMPLATE` branch commit history for a commit message like 'Merge branch 'dev' into TEMPLATE`'.

If you find evidence of this merge, rebuild your `TEMPLATE` branch from scratch:

## Rebuild the TEMPLATE branch

1. Clone the main nf-core pipeline repository to your local machine (not your development fork):

   ```bash
   git clone https://github.com/nf-core/<pipeline_name>.git
   cd <pipeline_name>
   ```

1. Retrieve the commit hash from when the original nf-core template was used to generate the pipeline (that is, with `nf-core pipelines create`). If you started with the nf-core template, check your git log:

   ```bash
   git checkout TEMPLATE
   git log --reverse
   ```

   The first commit typically represents the original template, with a commit message like `initial template build from nf-core/tools, version 1.9`.

1. Reset the `TEMPLATE` branch to this commit:

   ```bash
   git reset --hard <commit_hash>
   ```

1. Push the cleaned branch back to the repository:

   ```bash
   git push origin TEMPLATE --force
   ```

   :::warning
   This irreversibly overwrites git history.
Verify you have the correct branch name before running this command.
   :::

1. Switch back to `dev` and run the sync command to get the latest template version:

   ```bash
   git checkout dev
   nf-core pipelines sync
   ```

1. Delete your local copy of the pipeline that you cloned from the main nf-core repository.

1. Pull the fresh template branch into your personal fork:

   ```bash
   cd <path_to_forked_pipeline>
   git branch -D TEMPLATE
   git remote add upstream git@github.com:nf-core/<pipeline_name>.git
   git fetch upstream TEMPLATE
   git checkout --track upstream/TEMPLATE
   git push --force
   ```

You can now recreate the pull request from `TEMPLATE` into `dev`.
If conflicts occur, resolve them locally following the instructions in the sync documentation.

<!-- TODO: Add links to sync documentation -->
