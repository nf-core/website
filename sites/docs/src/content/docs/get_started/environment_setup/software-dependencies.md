---
title: Manage software dependencies
subtitle: Learn how to manage software dependencies
shortTitle: Software dependencies
weight: 3
---

Nextflow pipelines can specify software dependencies within the pipeline definition itself. Users are not required to manually install analysis tools, as the pipeline manages software provisioning through containerization or package management systems.

## How pipeline software works

nf-core pipelines are designed to work with multiple container and compute environment management systems. Each pipeline declares its required software tools, and users specify which management system to use via a configuration profile at runtime. The pipeline then automatically retrieves and utilizes the appropriate containers or environments for each step.

:::tip
The recommended software management system depends on your computing environment:

- **Local computers**: [Docker](#docker) is recommended for personal workstations and laptops.
- **HPC clusters**: Consult your system administrator, as many clusters have [Singularity](#singularity), [Apptainer](#apptainer), or [Shifter](#shifter) pre-installed.
- **When containers are unavailable**: [Conda or Mamba](#condamamba) provide fallback options.
  :::

## Software management systems

nf-core pipelines support several different systems for managing software dependencies. The best choice depends on your computing environment and whether you have administrator privileges.

### Docker

Docker is a containerization platform that packages software with all its dependencies into isolated containers. This platform is commonly deployed on local machines, single-user servers, and cloud computing platforms.

For installation instructions, see [Get Docker](https://docs.docker.com/get-docker/).

### Singularity

Singularity (now maintained as SingularityCE by Sylabs) is a container platform designed for high-performance computing (HPC) environments. Unlike Docker, it does not require administrator privileges to run containers, making it well-suited for shared computing clusters.

For installation instructions, see SingularityCE [Quick Start](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html).

### Apptainer

Apptainer is the open-source successor to Singularity, developed by the Linux Foundation after the project split in 2021. It provides similar functionality to Singularity and is becoming the standard on many HPC systems.

:::note
Apptainer currently uses Docker containers from nf-core rather than pre-built Apptainer/Singularity images.
:::

For installation instructions, see Apptainer [Quick start](https://apptainer.org/docs/user/latest/quick_start.html)

### Podman

Podman is a daemonless container engine that provides a Docker-compatible interface. It is designed as a drop-in replacement for Docker that does not require root privileges or a running daemon.

For installation instructions, see Podman [Installation instructions](https://podman.io/getting-started/installation)

### Charliecloud

Charliecloud is a lightweight container runtime designed specifically for HPC environments. It emphasizes simplicity and integration with existing HPC workflows.

For installation instructions, see Charliecloud [Installing](https://charliecloud.io/latest/install.html)

### Shifter

Shifter is a container runtime developed for HPC systems that allows Docker images to run on supercomputers and clusters. It is commonly deployed on large HPC facilities.

For installation instructions, see [Shifter](https://github.com/NERSC/shifter)

### Conda/Mamba

Conda and Mamba are package management systems that create isolated software environments. While they do not provide the same level of isolation as containers, they offer alternatives when containers are not available or practical.

:::caution
Conda environments may provide poorer reproducibility than container-based methods because low-level dependencies can change over time. Containers are recommended when possible.
:::

For installation instructions, see [Mamba installation](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)

<!-- TODO: Add additional resources -->
