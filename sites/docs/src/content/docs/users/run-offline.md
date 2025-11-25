---
title: Run pipelines offline
subtitle: Learn how to use nf-core pipelines without an internet connection
weight: 1
---

nf-core pipelines automatically fetche nearly everything it needs to run. However, you can also run nf-core pipelines on systems without internet access by preparing everything in advance.

You need:

- [Nextflow](#nextflow)
- [Pipeline code with its dependencies](#pipeline-code)
- [Reference genomes](#reference-genomes) (if applicable)

This guide walks you through running pipelines offline.

## Nextflow

To transfer Nextflow to an offline system:

1. [Install Nextflow](https://nextflow.io/docs/latest/getstarted.html#installation) in an online environment.
1. Run your pipeline locally.

    :::note
    Nextflow fetches the required plugins. It does not need to run to completion.
    :::

1. Copy the Nextflow binary and `$HOME/.nextflow` folder to your offline environment.
1. In your Nextflow configuration file, specify each plugin (both name and version), including default plugins.

    :::note
    This prevents Nextflow from trying to download newer versions of plugins.
    :::

1. Add the following environment variable in your `~/.bashrc` file:

   ```bash title=".bashrc"
   export NXF_OFFLINE='true'
   ```

## Pipeline code

To transfer pipeline code to an offline system:

1. Run `nf-core pipelines download <pipeline>` in an online environment.

   :::note
   This command will download the pipeline and configuration profiles, and packages everything into a `.tar.gz` file.
   :::

   :::tip
   Add the argument `--container singularity` to fetch the singularity container(s).
   :::

1. Transfer the `.tar.gz` file to your offline system and unpack it.

   :::note
   The archive contains directories called:

   - `workflow`: The pipeline files
   - `config`: [nf-core/configs](https://github.com/nf-core/configs) files
   - `singularity`: Singularity images (if you used `--container singularity`)

   The download tool adjusts the pipeline code to expect these relative paths. If you keep them together it should work as is.
   :::

   :::tip
   If you are downloading _directly_ to the offline storage (e.g., a head node with internet access whilst compute nodes are offline), use the `--singularity-cache-only` option for `nf-core pipelines download` and set the `$NXF_SINGULARITY_CACHEDIR` environment variable. This reduces total disk space by downloading singularity images to the `$NXF_SINGULARITY_CACHEDIR` folder without copying them into the target downloaded pipeline folder.
   :::

## Reference genomes

To use these references, download and transfer them to your offline cluster.
See [Reference genomes](./reference_genomes.md) to configure the references base path.

## Additional resources

For more information about running nf-core pipelines offline, see:

<!-- TODO: Add link to nf-core pipelines download documentation once available -->

- nf-core/bytesize: Running pipelines offline

  <!-- markdownlint-disable -->
  <iframe width="560" height="315" src="https://www.youtube.com/embed/N1rRr4J0Lps" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
  <!-- markdownlint-restore -->
