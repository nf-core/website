---
title: "Chapter 2: Getting started"
subtitle: "What you need before starting"
shortTitle: "Chapter 2: Getting started"
---

This chapter lists what you need before starting the training.

## Prerequisite knowledge

You should already be familiar with:

- Writing Nextflow pipelines (see the [Nextflow training](https://training.nextflow.io/)).
- How Nextflow modules and processes work (see the [Nextflow modules documentation](https://nextflow.io/docs/latest/module.html)).
- Git and GitHub basics (see the [git tutorial](https://git-scm.com/docs/gittutorial) and [GitHub Skills](https://skills.github.com/)).

## Required software and accounts

You need the following installed or set up:

- [Nextflow](https://www.nextflow.io/docs/latest/install.html).
- [nf-core tools](https://nf-co.re/docs/nf-core-tools/installation).
- A Nextflow-supported software management system, one of:
  - [Conda](https://conda-forge.org/download/)
  - [Docker](https://docs.docker.com/manuals/)
  - [Apptainer](https://apptainer.org/docs/user/latest/quick_start.html#installation)
- A [GitHub account](https://github.com).
- A fork and local clone of the [nf-core/modules repository](https://github.com/nf-core/modules).

## Tool assumptions

This training assumes you have:

- A command-line tool in mind to wrap as an nf-core module.
- A [Bioconda](https://bioconda.github.io/conda-package_index.html) recipe for that tool.

:::note
This training was written and tested with nf-core/tools 3.2.0 and Nextflow 24.10.4.
:::

The [next chapter](./3-what-is-a-module) introduces what defines an nf-core module.
