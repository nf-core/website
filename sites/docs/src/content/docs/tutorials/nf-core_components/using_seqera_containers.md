---
title: "Using Seqera Containers"
subtitle: "Step-by-step guide to using Seqera Containers in nf-core modules"
weight: 30
---

[Seqera Containers](https://seqera.io/containers) simplifies the container creation process for nf-core modules.
Currently, nf-core supports Seqera Containers as an alternative to the default biocontainers source. In the future, Seqera Containers will become the default source for all nf-core module containers.
For more information, see the [nf-core blog post](https://nf-co.re/blog/2024/seqera-containers-part-1).

Follow these steps to obtain both Docker and Singularity containers from the Seqera Containers website and configure your modules for integration with nf-core workflows.

0. (Optional but recommended) Test the `environment.yml` file locally to ensure that all packages and their versions are compatible:

   ```bash
   conda env create -f environment.yml
   ```

1. Open the [Seqera Containers website](https://seqera.io/containers/).
2. In the search bar, add all the required packages. You can copy and paste entries directly from your `environment.yml` file (for example, `conda-forge::<TOOL>=<VERSION>`) to ensure that the correct channels and versions are used.

   :::warning
   nf-core modules SHOULD include only a single tool and its required dependencies.
   Only add additional tools when necessary.
   :::

3. Set **Container setting** to `Docker`.
4. Select **Get Container**.
5. Copy the generated container URL and add it to the last line of the container definition in your module’s `main.nf` file.

   :::note
   The URL does not include `https://` at the beginning. Copy it exactly as shown.
   :::

6. Change the **Container settings** to `Singularity`.
7. Select **Get container**.
8. :grey_exclamation: Select the **HTTPS** checkbox.

   :::note
   This checkbox appears only after the container has been successfully built.
   :::

9. Copy the generated URL (including the `https://community-cr-prod...` prefix) and add it to the middle line of the container definition in your module’s `main.nf` file.
