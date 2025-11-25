---
title: Running nf-core pipelines offline
subtitle: Learn how to use nf-core pipelines without an internet connection
shortTitle: Running pipelines offline
weight: 1
---

nf-core pipelines are designed to automatically fetch everything they need to run, including pipeline code, container images, and reference data. However, many high-performance computing (HPC) environments restrict internet access on compute nodes for security reasons. In such cases, you can still run nf-core pipelines by preparing all required components in advance on a system with internet access, then transferring them to your offline environment.

Running pipelines offline requires three main components:

- **Nextflow runtime**: The Nextflow workflow engine and its plugins
- **Pipeline code**: The nf-core pipeline itself, along with container images
- **Reference data**: Genome files and other reference datasets (if needed by your pipeline)

## Offline setup

The following sections describe how to prepare each component for offline execution.

### Transfer Nextflow offline

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

### Transfer pipeline code offline

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
     :::

   :::tip
   If you are downloading _directly_ to the offline storage (e.g., a head node with internet access whilst compute nodes are offline), use the `--singularity-cache-only` option for `nf-core pipelines download` and set the `$NXF_SINGULARITY_CACHEDIR` environment variable. This reduces total disk space by downloading singularity images to the `$NXF_SINGULARITY_CACHEDIR` folder without copying them into the target downloaded pipeline folder.
   :::

### Transfer reference genomes offline

To use nf-core reference genomes offline, download and transfer them to your offline cluster.
See [Reference genomes](./reference_genomes.md) for more information.

## Additional resources

For more information about running nf-core pipelines offline, see:

<!-- TODO: Add link to nf-core pipelines download documentation once available -->

- nf-core/bytesize: Running pipelines offline

  <!-- markdownlint-disable -->
  <iframe width="560" height="315" src="https://www.youtube.com/embed/N1rRr4J0Lps" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
  <!-- markdownlint-restore -->
