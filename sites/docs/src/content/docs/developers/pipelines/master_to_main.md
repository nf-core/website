---
title: Renaming branches
subtitle: How to switch your pipeline's default branch from `master` to `main`
shortTitle: Renaming branches
---

Pipeline developers may choose to change the pipeline's default branch from `master` to `main` for consistency with other repositories.

This guide demonstrates how to change a pipeline's default branch from `master` to `main`.

:::note{title="Prerequisites"}
You will need the following to get started:

- Admin privileges for the GitHub repository
  - If you don’t have permission, ask `@core-team` for help on Slack.
- nf-core tools

:::

## Rename `master` to `main`

1. Open your repository on GitHub while viewing the `master` branch
1. Select **Branches** (located to the right of the branch dropdown menu)
1. In the **Default** branch section, find `master`, select the three-dot menu, then choose **Rename branch**
1. Rename `master` to `main`
1. Select **learn more** and save the displayed instructions for contributors
1. Select **Rename branch** to confirm
1. Verify you're now on `main` in the **Code** tab
1. In your IDE, switch to the `dev` branch
1. Run the following command to obtain the `main` branch locally:

   ```bash
   git fetch upstream
   ```

1. Configure Git to use `main` as the default branch name:

   ```bash
   git config --global init.defaultBranch main
   ```

1. Verify the configuration:

   ```bash
   git config --global init.defaultBranch
   ```

1. While on `dev`, run the sync command:

   ```bash
   nf-core pipelines sync
   ```

1. Create a new branch for your changes:

   ```bash
   git switch -c default-branch-change
   ```

1. Merge the template changes:

   ```bash
   git merge TEMPLATE
   ```

   :::note
   Ensure your `nf-core/tools` version matches the template version in your pipeline.
   :::

1. Resolve any merge conflicts that appear

   :::tip
   If it’s the ROcrate JSON file, you can accept all incoming change.
   :::

1. Search your repository for references to `master` and update them to `main`

   :::warning
   Make sure not to modify references of master in links to other repositories! If in doubt, ask on the nf-core Slack!
   :::

1. Run the linting tool to check for issues:

   ```bash
   nf-core pipelines lint
   ```

1. Commit your changes:

   ```bash
   git add -am 'Change default branch'
   ```

1. Push the changes to GitHub
1. Create a pull request against `dev` on GitHub
1. Request community review on Slack
1. Merge the pull request upon approval

## Update collaborator repositories

After merging the changes, notify your collaborators to update their clones and forks. They should run the following commands:

```bash
git branch -m master main
git fetch origin
git branch -u origin/main main
git remote set-head origin -a
```

You should also update your own clones and forks using these commands.
