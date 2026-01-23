---
title: Template syncs
subtitle: How nf-core pipelines are kept up to date with community standards
weight: 1
---

nf-core maintains community standards across all pipelines through automated template synchronisation.
When the nf-core template updates with new features or best practices, these changes automatically propagate to every pipeline in the community.
This automation reduces maintenance burden for pipeline developers and ensures consistent quality across the nf-core ecosystem.

## How template syncs work

Template syncs use git to track changes between your pipeline and the base template.
Each pipeline repository includes a special `TEMPLATE` branch that stores the unmodified template code.
When nf-core releases a new tools version, GitHub Actions automatically creates a pull request for each pipeline.
These PRs show only the template changes relevant to your pipeline, excluding your custom modifications.
You review and merge these PRs to keep your pipeline aligned with community standards.

The sync requires that your `TEMPLATE` branch shares git history with your `dev` branch.
The `nf-core pipelines create` command sets this up automatically.
If you created your pipeline another way, you need to configure this manually.

## Working with template syncs

This section covers common template sync scenarios and how to handle them.

### Merge automated PRs

When nf-core releases a new tools version, automated pull requests are created for each pipeline.
Learn how to review and merge these PRs, whether they contain simple changes or complex merge conflicts.
This guide covers both web-based conflict resolution for minor conflicts and local resolution for major conflicts.

See [Merge automated PRs](merge-automated-pull-requests.md) for more information.

### Manually sync your pipeline

In rare cases, you need to trigger synchronisation manually.
This applies when automated sync did not run or when you maintain a custom pipeline.
Learn how to run the sync command manually, either for official nf-core pipelines or custom pipelines with pull request creation.

See [Manually sync your pipeline](manual-sync.md) for more information.

### Fix a broken TEMPLATE branch

If you resolve merge conflicts through the GitHub web interface, your `TEMPLATE` branch can become corrupted with `dev` branch history.
This guide shows you how to identify a broken `TEMPLATE` branch and rebuild it from scratch to restore proper template syncing.

See [Fix a broken TEMPLATE branch](fix-broken-template-branch.md) for more information.

### Set up a pipeline sync retrospectively

If your pipeline was not created with `nf-core pipelines create`, your `TEMPLATE` branch may not be configured correctly.
This guide walks you through creating a proper `TEMPLATE` branch and merging it into your existing pipeline to enable future template syncs.

See [Set up a pipeline sync retrospectively](set-up-pipeline-sync.md) for more information.
