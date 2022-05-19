---
title: 3 - Configuration
---

## Configuration of Gitpod

Within each git repository, the main file that controls the gitpod environment is the `.gitpod.yml` file, which contains the instructions on which environment to build and which tools to install.

For example, check out the nf-core modules `.gitpod.yml` file [here](https://github.com/nf-core/modules/blob/master/.gitpod.yml). You majority of nf-core gitpod ymls look like this, with two sections:

1. **image**: which pulls the nfcore image `nfcore/gitpod:latest` into Gitpod., which contains `nextflow` and nf-core scripts and other essential tools such as docker.

2. **vscode**: this allows vscode extensions within your environment.

## Configuring your own custom gitpod environment

There are a few more common adaptations to `./gitpod.yml` files, which include:

1. **github** - allows configuration of github. e.g. Allows gitpod to create prebuilds for branches.
2. **ports** - opens a port to serve traffic to a public URL. This is used within the nfcore webpage yml (see: https://github.com/nf-core/nf-co.re/blob/master/.gitpod.yml).
3. **tasks** - this tells gitpod to run particular jobs. Usually you will see the following:

`- init:` sections can be used to install packages as a pre-build, so it doesn't have to run each time you open an environment.

`- command:` sections execute the given lines of code on every workspace startup.

For more detailed information about these settings, check out the extensive docs at Gitpod [here](https://www.gitpod.io/docs/config-gitpod-file).
