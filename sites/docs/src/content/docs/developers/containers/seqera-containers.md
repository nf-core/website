---
title: Seqera Containers
subtitle: Using Seqera Containers in nf-core modules
shortTitle: Seqera Containers
weight: 2
---

Seqera Containers simplifies the container creation process for nf-core modules. This guide explains how to obtain Docker and Singularity containers from Seqera Containers for your nf-core modules.

:::note
Currently, nf-core supports Seqera Containers as an alternative to the default biocontainers source. In the future, Seqera Containers will become the default source for all nf-core module containers. For more information, see the [nf-core blog post](https://nf-co.re/blog).
:::

:::note{title="Prerequisites"}
You will need the following to get started:

- An `environment.yml` file defining your module's dependencies
- The package names for your required tools
- Access to [Seqera Containers](https://seqera.io/containers/)

:::

## Module package requirements

nf-core modules should include only a single tool and its required dependencies. Avoid bundling multiple unrelated tools in a single container.

## Generate containers

Follow these steps to generate both Docker and Singularity container URLs for your module:

1. Visit the [Seqera Containers website](https://seqera.io/containers/)

1. Add your required packages to the search bar. You can paste package names directly from your `environment.yml` file.

    :::tip
    Test your `environment.yml` file locally using conda before generating containers. This helps verify package compatibility and catch potential issues early.
    :::

1. Set the **Container** setting to **Docker**

1. Select **Get Container**

1. Copy the generated Docker URL (without the `https://` prefix)

1. Change the **Container** setting to **Singularity**

1. Select **Get container**

1. Wait for the build to complete successfully

1. Check the **HTTPS** checkbox

    - This option only appears after a successful build

1. Copy the generated Singularity URL including the full `https://community-cr-prod...` prefix

1. Add both URLs to your module's container definition in `main.nf`:

   ```groovy
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'https://community-cr-prod.seqera.io/<singularity-container-id>' :
       'community.wave.seqera.io/<docker-container-id>' }"
   ```

    :::warning
    Docker URLs should not include the `https://` prefix, while Singularity URLs must include it. The HTTPS checkbox for Singularity only appears after the container builds successfully.
    :::

## Additional resources

- [Seqera Containers website](https://seqera.io/containers/)
- [Contributing components guide](../../contributors/contribute-components.md)
- [nf-core/modules repository](https://github.com/nf-core/modules)
