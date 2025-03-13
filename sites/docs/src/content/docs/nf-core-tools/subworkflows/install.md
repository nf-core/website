---
title: Install a subworkflow
subtitle: Install subworkflows in a pipeline
shortTitle: install
weight: 30
---

You can install subworkflows from [nf-core/modules](https://github.com/nf-core/modules) in your pipeline using `nf-core subworkflows install`.
A subworkflow installed this way will be installed to the `./subworkflows/nf-core` directory.

<!-- RICH-CODEX
working_dir: tmp/nf-core-nextbigthing
before_command: >
  echo "repository_type: pipeline" >> .nf-core.yml
-->

![`nf-core subworkflows install bam_rseqc`](/images/tools/nf-core-subworkflows-install.svg)

You can pass the subworkflow name as an optional argument to `nf-core subworkflows install` like above or select it from a list of available subworkflows by only running `nf-core subworkflows install`.

There are four additional flags that you can use when installing a subworkflow:

- `--dir`: Pipeline directory, the default is the current working directory.
- `--force`: Overwrite a previously installed version of the subworkflow.
- `--prompt`: Select the subworkflow version using a cli prompt.
- `--sha <commit_sha>`: Install the subworkflow at a specific commit.
