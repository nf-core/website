---
title: Environment setup
subtitle: Setting up your environment for nf-core
parentWeight: 2
weight: 1
---

To work with nf-core pipelines, you will need to set up your computing environment with some essential tools. This section guides you through installing and configuring everything you need to get started.

## What you'll need

<!-- TODO: Fix astro errors with links (Cannot read properties of undefined (reading 'id')) -->

The core requirements for running nf-core pipelines are:

1. **[Nextflow](./nextflow.md)**: The workflow management system that powers all nf-core pipelines
2. **[Container or compute environment manager](./software-dependencies.md)**: A system to handle pipeline software (e.g., Docker, Singularity, or Conda)

These additional tools are recommended for development:

3. **[nf-core tools](./nf-core-tools.md)**: Command-line tools for working with nf-core pipelines (optional, but recommended)
4. **[VS Code](./vs-code.md)**: A code editor with Nextflow extensions (optional, for pipeline development)

## Getting started

If you are new to nf-core, we recommend following the guides in order:

1. [Install Nextflow](./nextflow.md)
2. Set up a [Container or compute environment manager](./software-dependencies.md) based on your system
3. Install [nf-core tools](./nf-core-tools.md) for additional utilities and pipeline management features
4. Set up [VS Code](./vs-code.md) with the recommended extensions ff you plan to develop pipelines

:::note
System requirements and installation methods vary depending on your operating system and whether you have administrator privileges. If you are uncertain, speak to your system administrator.
:::
