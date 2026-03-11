---
title: "nf-core/tools - 3.4.0"
subtitle: "Container regex-it"
headerImage: "https://images.unsplash.com/photo-1592963219683-3593d83276af"
headerImageAlt: "Photo of a man staring at a high stack of shipping containers."
pubDate: 2025-10-16T12:00:00+02:00
authors:
  - "mashehu"
  - "mirpedrol"
label:
  - "tools"
maxHeadingDepth: 3
---

This release brings significant improvements to the download command, a new devcontainer setup, and better ARM64 architecture handling.
As always, if you have any problems or run into any bugs, reach out on the [#tools slack channel](https://nfcore.slack.com/archives/CE5LG7WMB).

# Highlights

- [Refactored download command](#refactored-download-command)
- [New devcontainer setup](#new-devcontainer-setup)
- [Improvements in ARM64 architecture handling](#improvements-in-arm64-architecture-handling)
- [CLI convenience improvements](#cli-convenience-improvements)

## Refactored download command

The `nf-core download` command has received a major overhaul (see the [blog post](/blog/2025/refurbushing-the-pipeline-download) for more details about the motivation and new approach), bringing several powerful new features and improvements:

#### Nextflow inspect for container discovery

The download command now uses `nextflow inspect` to discover containers used in pipelines, replacing the legacy regex-based approach.
This provides more accurate and reliable container detection.

:::note
This requires Nextflow version 25.04.04 or later.
:::

The old regex-based container discovery could sometimes miss containers or pick up false positives.
With `nextflow inspect`, Nextflow itself tells us exactly which containers are used, making the process more robust and future-proof.

#### Docker tar archive support

You can now download Docker images directly into tar archives, making it easier to transfer and deploy pipelines in air-gapped environments:

```bash /--compress tar/
nf-core download <pipeline> --container-system docker --compress tar
```

## New devcontainer setup

With the [sunsetting of Gitpod](https://ona.com/stories/gitpod-classic-payg-sunset) (😢) we recommend using GitHub Codespaces instead.

We updated therefore the devcontainer configuration to make the experience as seamless as possible.

:::note
If you want to run a Nextflow pipeline, we recommend to use `-profile singularity` when running in a codespace. Docker is not fully supported yet, and it can lead to unexpected behaviours.
:::

## Improvements in ARM64 architecture handling

The template now includes better support for ARM64 architectures:

- The generic `arm` profile has been replaced with a more specific `arm64` profile
- A new `emulate_amd64` profile has been added for running AMD64 containers on ARM64 systems using emulation, which is better suited for Apple Silicon.

This provides clearer options for pipeline users working with Apple Silicon Macs and ARM64 servers

## CLI convenience improvements

### `modules lint` and `modules bump-versions` all sub-tools

If you want to bump the version of a series of modules which are based on the same tool, e.g., samtools, you can now just specify the tool name and the same command will be run on all sub-tools:

```bash
nf-core modules lint samtools
nf-core bump-versions samtools
```

Thank you [@nh13](https://github.com/nh13) for adding this feature!

### Command short-hands

All main commands have now shorter aliases:

- `nf-core pipelines{:bash}` -> `nf-core p{:bash}`
- `nf-core modules{:bash}` -> `nf-core m{:bash}`
- `nf-core subworkflows{:bash}` -> `nf-core s{:bash}`
- `test-datasets{:bash}` -> `nf-core tds{:bash}`

Thanks to the new version of [rich-click](https://ewels.github.io/rich-click/1.9/blog/2025/09/16/version-1.9/) for adding this feature (and nicer themeing)!

# Pipeline template changes

This release includes also some changes to the pipeline template.

### CI

- We switched the `download_pipeline.yml` action to use the latest version of nf-core/tools instead of `dev`.
- We fixed the failing `aws{full}test.yml` action by switching to organization-wide variables as inputs.

### Release announcements

Mastodon announcements now include the pipeline description, providing more context when pipelines are released.

# Miscellaneous

## Modules linting

Modules can now use Nextflow's `exec:` blocks and we lint for possible conflicts with `shell:` blocks.

Thanks to [@muffato](https://github.com/muffato) for implementing this feature!

## Linting

The linter now uses the organization specified in `.nf-core.yml` when checking manifest names, homepages, and MultiQC config comments.
This makes it easier to fork and customize nf-core pipelines for your own organization while maintaining proper attribution.

Thank you [@rrahn](https://github.com/rrahn) for implementing this feature!

## Updated dependencies

- **python**: 3.9 reached end of life, so we set 3.10 as the minimum version and also expanded to support Python 3.14
- **nf-schema**: bumped to 2.5.0 with improved help message creation for future Nextflow versions
- **Minimum Nextflow version bumped to 25.04.0** to support the new download logic

# Changelog

You can find the complete changelog and technical details [on GitHub](https://github.com/nf-core/tools/releases/tag/3.4.0).

# Patch release 3.4.1

We found two small bugs (one whitespace error in `nextflow.config` and the devcontainer setup for pipelines was faulty) and fixed them for 3.4.1.

# Getting help to update your pipeline template

All the template changes are listed in [the Template section of the changelog](https://github.com/nf-core/tools/releases/tag/3.4.0).
Below, we have collected the most common merge conflicts that you may find and how to resolve them.

### \*.nf.test.snap - nf-test snapshots

##### Changes

If you use MultiQC in your pipeline, your nf-test snapshots will be outdated, because we updated the MultiQC version to 1.31.

##### Resolution

Regenerate the snapshots.

### nextflow.config

##### Changes

- The generic `arm` profile has been replaced with a more specific `arm64` and `emulate_amd64` profile for better ARM64 architecture handling.
- We also bumped the Nextflow version.
- the gitpod profile has been removed

##### Resolution

- If you had custom changes in your `arm` profile, move them over to the `arm64` profile, otherwise accept the changes from the template.
- Accept the new Nextflow version (if yours is not higher) and keep
- Accept the removed gitpod profile

### .github/workflows/awstest.yml and .github/workflows/awsfulltest.yml

##### Changes

We replaced many of the `secrets.` variables with organization-wide `vars.` variables.

##### Resolution

Accept the new variable references (`vars.` instead of `secrets.`).
If you previously changed one of these values or added custom changes to these workflows, keep both the new variable references and your modifications.

### .nftignore

##### Changes

Changes related to the MultiQC update have been made to this file.

##### Resolution

Accept all changes from the template and keep any additional patterns you added.
Be sure not to add an empty line otherwise pre-commit will complain

### main.nf

##### Changes

New inputs have been added and the `PIPELINE_INITIALIZATION` subworkflow has been updated.

##### Resolution

Accept new inputs and ensure the `PIPELINE_INITIALIZATION` subworkflow changes are integrated.
Keep any custom logic you added while merging the template updates.

### .gitpod.yml

##### Changes

Gitpod configuration has been removed in favor of GitHub Codespaces with devcontainers.

##### Resolution

You can safely delete the `.gitpod.yml` file.

### subworkflows/local/utils*nfcore*\*\_pipeline/main.nf

##### Changes

Utility subworkflows have been updated.

##### Resolution

Accept changes from the template for the local utility subworkflows to ensure compatibility with the latest nf-core standards.
