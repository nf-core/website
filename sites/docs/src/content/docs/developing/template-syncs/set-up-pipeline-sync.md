---
title: Setting up a pipeline sync retrospectively
subtitle: How to set up a correct TEMPLATE branch in your nf-core pipeline repository
shortTitle: Setting up a pipeline sync retrospectively
weight: 5
---

Use this guide to set up a correct `TEMPLATE` branch for pipelines that were not created with `nf-core pipelines create`.
If you created your pipeline with `nf-core pipelines create`, your `TEMPLATE` branch is already configured correctly.

:::warning
Before proceeding, push all local changes to GitHub and consider making a backup clone of your repository.
:::

Alternatively, consider restarting your pipeline project by running `nf-core pipelines create` and copying your modifications into the newly created pipeline.

:::note{title="Prerequisites"}
This guide assumes you are working directly with the main nf-core repository.
You can also work on your own fork.
:::

## Create the TEMPLATE branch

1. Clone your pipeline into a new directory:

   ```bash
   mkdir <temp_directory>
   cd <temp_directory>
   git clone https://github.com/nf-core/<pipeline_name>.git
   cd <pipeline_name>
   ```

1. Create the new `TEMPLATE` branch and delete all files to create a completely empty branch:

   ```bash
   git checkout --orphan TEMPLATE && git rm -rf '*'
   ```

1. Verify your branch is completely empty:

   ```bash
   git status
   ```

   You should see:

   ```console
   On branch TEMPLATE

   No commits yet

   nothing to commit (create/copy files and use "git add" to track)
   ```

1. Regenerate your pipeline from scratch using the most recent template:

   ```bash
   nf-core pipelines create --no-git
   ```

   If your pipeline already has versioned releases (for example, you are not currently on `1.0dev`), specify the version number:

   ```bash
   nf-core pipelines create --no-git --version 1.3dev
   ```

   :::note
   The version you choose should match the branch you intend to merge with.
If you already have a release, use the version number specified in your `dev` branch.
   :::

1. Follow the prompts to enter the pipeline name, description, and authors. Use the exact text from your existing `nextflow.config` file (`manifest.name` etc.).

1. Move the newly created template files into your root git directory:

   ```bash
   mv nf-core-<pipeline_name>/* .
   mv nf-core-<pipeline_name>/.[!.]* .
   rmdir nf-core-<pipeline_name>
   ```

1. Verify the newly created files are in the correct place:

   ```bash
   git status
   ```

   You should see output similar to:

   ```console
   On branch TEMPLATE

   No commits yet

   Untracked files:
     (use "git add <file>..." to include in what will be committed)

     .editorconfig
     .gitattributes
     .github/
     .gitignore
     .markdownlint.yml
     CHANGELOG.md
     CITATIONS.md
     CODE_OF_CONDUCT.md
     LICENSE
     README.md
     assets/
     bin/
     conf/
     docs/
     lib/
     main.nf
     modules.json
     modules/
     nextflow.config
     nextflow_schema.json
     subworkflows/
     workflows/

   nothing added to commit but untracked files present (use "git add" to track)
   ```

1. Commit the template files:

   ```bash
   git add .
   git commit -m "Initial template commit"
   ```

1. Push the `TEMPLATE` branch to the upstream nf-core repository:

   ```bash
   git push --set-upstream origin TEMPLATE
   ```

## Merge TEMPLATE into main branches

Merge the `TEMPLATE` branch into your main pipeline branches.
This requires manually resolving all merge conflicts.

1. Check out your development branch and create a new branch for the merge:

   ```bash
   git checkout dev
   git checkout -b template_merge
   ```

1. Merge the `TEMPLATE` branch:

   ```bash
   git merge TEMPLATE --allow-unrelated-histories
   ```

   :::note
   If the merge command shows entire files as new, try adding the `-Xignore-space-at-eol` flag.
   :::

1. Resolve merge conflicts. You may see many conflicts:

   ```console
   Auto-merging nextflow.config
   CONFLICT (add/add): Merge conflict in nextflow.config
   Auto-merging main.nf
   CONFLICT (add/add): Merge conflict in main.nf
   Auto-merging environment.yml
   CONFLICT (add/add): Merge conflict in environment.yml
   Auto-merging docs/usage.md
   CONFLICT (add/add): Merge conflict in docs/usage.md
   ```

   Go through each file to resolve the merge conflicts.
Use a visual merge tool to avoid mistakes when handling many merge markers.

1. Commit and push the resolved changes:

   ```bash
   git commit -m "Merged vanilla TEMPLATE branch into main pipeline"
   git push origin template_merge
   ```

1. Create a pull request from your `template_merge` branch to the main nf-core repository.

1. Once merged into the `dev` branch, future automatic template syncs will work correctly. When nf-core releases new versions of `nf-core/tools`, pull requests will automatically be created to merge updates into your pipeline.
