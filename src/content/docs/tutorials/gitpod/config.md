---
title: Configuration
subtitle: Configuring your Gitpod environment
weight: 3
---

## Configuration of Gitpod

Within each git repository, the main file that controls the Gitpod environment is the `.gitpod.yml` file, which contains the instructions on which environment to build and which tools to install.

For example, check out the nf-core modules `.gitpod.yml` file [here](https://github.com/nf-core/modules/blob/master/.gitpod.yml). The majority of nf-core Gitpod YAML files look like this, with two sections:

1. **image**: this pulls the nf-core image `nfcore/gitpod:latest` into Gitpod, which contains `nextflow` and nf-core scripts and other essential tools such as Docker.

2. **vscode**: this allows the source-code editor vscode extensions within your environment.

## Configuring your own custom Gitpod environment

There are a few more common adaptations to `./gitpod.yml` files, which include:

1. **github**: this allows configuration of Github. e.g. allowing Gitpod to create prebuilds for branches.
2. **ports**: this opens a port to serve traffic to a public URL. This is used within the nf-core webpage yml (see: https://github.com/nf-core/website/blob/main/.gitpod.yml).
3. **tasks**: this tells Gitpod to run particular jobs. Usually you will see the following:

`- init:` sections that can be used to install packages as a pre-build, so it doesn't have to run each time you open an environment.

`- command:` sections that execute the given lines of code on every workspace startup.

For more detailed information about these settings, check out the extensive docs at Gitpod [here](https://www.gitpod.io/docs/config-gitpod-file).
