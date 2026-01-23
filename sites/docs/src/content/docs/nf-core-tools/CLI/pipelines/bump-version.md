---
title: Bump version number
subtitle: Bump a pipeline version number
shortTitle: bump-version
weight: 90
---

When you release a new version of a nf-core pipeline, you need to update version numbers in several places.
The `nf-core pipelines bump-version` command automates this process to avoid manual errors.

This command requires your pipeline to pass the linting tests.

Run `nf-core pipelines bump-version <new_version>`, for example:

<!-- RICH-CODEX
working_dir: tmp/nf-core-nextbigthing
-->

![`nf-core pipelines bump-version 1.1`](../../../../../assets/images/tools/nf-core-bump-version.svg)

To specify a different directory, use `--dir <pipeline_dir>`.
To change the required version of Nextflow instead of the pipeline version number, use the `--nextflow` flag.
