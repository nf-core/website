---
title: Bumping version number
subtitle: Bumping a pipeline version number
shortTitle: bump-version
weight: 90
---

When releasing a new version of a nf-core pipeline, version numbers have to be updated in several different places. The helper command `nf-core pipelines bump-version` automates this for you to avoid manual errors (and frustration!).

The command uses results from the linting process, so will only work with workflows that pass these tests.

Usage is `nf-core pipelines bump-version <new_version>`, eg:

<!-- RICH-CODEX
working_dir: tmp/nf-core-nextbigthing
-->

![`nf-core pipelines bump-version 1.1`](../../../../assets/images/tools/nf-core-bump-version.svg)

You can change the directory from the current working directory by specifying `--dir <pipeline_dir>`. To change the required version of Nextflow instead of the pipeline version number, use the flag `--nextflow`.
