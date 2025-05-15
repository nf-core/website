---
title: "Using Seqera Containers"
subtitle: "Step by step guide to using Seqera Containers in nf-core modules"
weight: 30
---
[Seqera Containers](https://seqera.io/containers) simplifies the container creation process for nf-core modules. The following steps will guide you through obtaining both Docker and Singularity containers from the Seqera website to properly configure your module for integration with nf-core workflows.
0. (Optional, but recommended) test building the `environment.yml` locally to ensure that the packages and their versions are compatible:

    ```bash
    conda env create -f environment.yml
    ```

1. Go to [https://seqera.io/containers/](https://seqera.io/containers/)
2. Add all the packages in the search bar (you can copy paste e.g. `conda-forge::<tool>=<version>` from your `environment.yml` file, to make sure you have exactly the right channels/versions).
3. Make sure container setting is on `Docker`.
4. Press 'Get Container'.
5. Copy the displayed URL and put it on the last line of your module's `main.nf` container definition (this does _not_ have a `https://` at the beginning!).
6. Repeat steps 2-5.
7. Switch the 'Container settings' to `Singularity`.
8. Press 'Get container'.
9. :grey_exclamation: Tick the `HTTPS` checkbox.
10. Copy the URL (_with_ the `https://community-cr-prod<...>` at the beginning), and put it on the middle line of your module's `main.nf` definition.
11. Repeat steps 8-10.

:::note
See the blog post about our [Migration from Biocontainers to Seqera Containers](https://nf-co.re/blog/2024/seqera-containers-part-1) for more information on why we are using Seqera Containers and how they compare to Biocontainers.
:::
