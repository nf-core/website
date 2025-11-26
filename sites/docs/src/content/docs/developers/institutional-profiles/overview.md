---
title: Institutional profiles
subtitle: Learn how to write and contribute institutional configuration profiles
weight: 1
---

nf-core offers centralized storage for Nextflow configuration profiles at the institutional level. These profiles enable sharing of common pipeline configurations across all users of a particular cluster and can be used in any Nextflow pipeline, not just nf-core pipelines.

## What are institutional profiles?

Institutional profiles are configuration files that define cluster-specific settings for running Nextflow pipelines. They specify details such as:

- Job scheduler settings (SLURM, SGE, PBS, etc.)
- Resource limits (maximum CPUs, memory, walltime)
- Queue or partition configurations
- Container engine settings (Singularity, Docker, Podman, Conda)
- Shared resource paths (reference genomes, container caches)

By centralizing these configurations in the [nf-core/configs](https://github.com/nf-core/configs) repository, all users at an institution can benefit from tested, optimized settings without maintaining separate configuration files.

:::note
While institutional profiles are best contributed to [nf-core/configs](https://github.com/nf-core/configs) for community benefit, this documentation can also be used to write well-structured configuration profiles for private or personal use.
:::

## Profile types

Two types of institutional profiles exist:

- **Global profiles**: Apply to all nf-core pipelines, defining cluster-level settings that work across all workflows
- **Pipeline-specific profiles**: Apply to individual nf-core pipelines with customized parameters for specific workflows

Most institutions only need global profiles. Pipeline-specific profiles are useful when certain pipelines require unique settings due to unusual resource requirements or dependencies.

## How users access profiles

Once your profile is merged into the nf-core/configs repository, users can access it by specifying the profile name with the `-profile` flag:

```bash
nextflow run nf-core/rnaseq -profile <your_institution>
```

Users can combine institutional profiles with other profiles using comma-separated values:

```bash
nextflow run nf-core/rnaseq -profile <your_institution>,docker
```

## Writing an institutional profile

This section guides you through creating and contributing an institutional profile to the nf-core community.

### Gather cluster information

Before writing your profile, collect essential information about your cluster's specifications, scheduler, resource limits, and available tools. Understanding these details ensures your profile works correctly for all users.

See [Preparing to write a profile](./preparing-to-write.md) for more information.

### Understand file structure

Institutional profiles require specific files in the nf-core/configs repository. Learn which files you need to create or modify to add your profile.

See [File structure requirements](./file-structure.md) for more information.

### Configure profile settings

Configuration files use Nextflow's configuration syntax to define parameters, process settings, executor options, and container configurations. Understanding each configuration scope helps you create effective profiles.

See [Configuration file components](./configuration.md) for more information.

### Document your profile

Each profile requires documentation that explains its purpose, usage, and any special requirements. Clear documentation helps users at your institution get started quickly.

See [Documentation requirements](./documentation.md) for more information.

### Test your profile

Before submitting your profile, test it thoroughly to verify correct settings, container engine detection, and scheduler job submissions. Testing prevents issues for future users.

See [Testing profiles](./testing.md) for more information.
