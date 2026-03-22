---
title: Environment setup
subtitle: Set up your environment
parentWeight: 2
weight: 1
---

To work with nf-core pipelines, you will need to set up your computing environment with some essential tools.
This section guides you through installing and configuring everything you need to get started.

## What you'll need

<!-- TODO: Fix astro errors with links (Cannot read properties of undefined (reading id)) -->

The core requirements for running nf-core pipelines are:

1. **[Nextflow](./nextflow)**: The workflow management system that powers all nf-core pipelines
1. **[Container or compute environment manager](./software-dependencies)**: A system to handle pipeline software (e.g., Docker, Singularity, or Conda)

These additional tools are recommended for development:

1. **[nf-core tools](./nf-core-tools)**: Command-line tools for working with nf-core pipelines (optional, but recommended)
1. **[VS Code](./vs-code)**: A code editor with an official Nextflow extension (optional, for pipeline development)
1. **[Prettier](./prettier)**: A code styling autoformatter used in multiple nf-core repositories (optional, but recommended)

## Getting started

If you are new to nf-core, we recommend following the guides in order:

1. [Install Nextflow](../environment_setup/nextflow)
1. Set up a [Container or compute environment manager](../environment_setup/software-dependencies) based on your system

If you then decide to develop with the community

1. Install [nf-core tools](../environment_setup/nf-core-tools) for additional utilities and pipeline management features
1. Set up [VS Code](../environment_setup/vs-code) with the recommended extensions if you plan to develop pipelines
1. Set up [Prettier](../environment_setup/prettier) so you can have auto-formatting of your code before pushing to GitHub

:::note
System requirements and installation methods vary depending on your operating system and whether you have administrator privileges. If you are uncertain, speak to your system administrator.
:::

:::tip
nf-core repositories include Development Container configurations that provide ready-to-use environments with all required software, accessible via GitHub Codespaces or VS Code.
See [Dev Containers](../environment_setup/dev-containers) for more information.
:::

## Need help?

If you encounter issues during setup, check our [troubleshooting guide](../environment_setup/troubleshooting) for solutions to common problems with Nextflow, nf-core tools, and development environments.

Still stuck? The nf-core community is here to help! Join us on [Slack](https://nf-co.re/join) where you can ask questions, get support, and connect with other users and developers.
