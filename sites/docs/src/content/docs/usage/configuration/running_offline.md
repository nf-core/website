---
title: Running offline
subtitle: Run nf-core pipelines offline
shortTitle: Running offline
weight: 4
---

## Running offline

When Nextflow is connected to the internet it will fetch nearly everything it needs to run a pipeline. Nextflow can also run analysis on an offline system that has no internet connection. However, there are a few extra steps that are required to get everything you will need locally.

To run a pipeline offline you will need three things:

- [Nextflow](#nextflow)
- [Pipeline assets](#pipeline-assets)
- [Reference genomes](#reference-genomes) _(if required)_

These will first need to be fetched on a system that _does_ have an internet connection and transferred to your offline system.

### Nextflow

You need to have Nextflow installed on your local system.
You can do this by installing Nextflow on a machine that _does_ have an internet connection and transferring to the offline system:

1. Install Nextflow locally
   :::warning
   do _not_ use the `-all` package, as this does not allow the use of custom plugins.
   :::
2. Run a Nextflow pipeline locally so that Nextflow fetches the required plugins.
   - It does not need to run to completion.
3. Copy the Nextflow executable and your `$HOME/.nextflow` folder to your offline environment
4. Specify the name and version each plugin that you downloaded in your local Nextflow configuration file
   - This will prevent Nextflow from trying to download newer versions of plugins.
5. Set `export NXF_OFFLINE='true'`
   - Add this command to your shell configuration file (e.g., `~/.bashrc` or `~/.zshrc`) to add this directory to `PATH` permanently

### Pipeline assets

To run a pipeline offline, you need the pipeline code, the software dependencies, and the shared nf-core/configs profiles.
We have created a helper tool as part of the _nf-core_ package to automate this for you.

On a computer with an internet connection, run `nf-core download <pipeline>` to download the pipeline and config profiles.
Add the argument `--container singularity` to also fetch the singularity container(s).

The pipeline and requirements will be downloaded, configured with their relative paths, and packaged into a `.tar.gz` file by default.
This can then be transferred to your offline system and unpacked.

Inside, you will see directories called `workflow` (the pipeline files), `config` (a copy of [nf-core/configs](https://github.com/nf-core/configs)), and (if you used `--container singularity`) a directory called `singularity`.
The pipeline code is adjusted by the download tool to expect these relative paths, so as long as you keep them together it should work as is.

### Shared storage

If you are downloading _directly_ to the offline storage (e.g., a head node with internet access whilst compute nodes are offline), you can use the `--singularity-cache-only` option for `nf-core download` and set the `$NXF_SINGULARITY_CACHEDIR` environment variable.
This downloads the singularity images to the `$NXF_SINGULARITY_CACHEDIR` folder and does not copy them into the target downloaded pipeline folder.
This reduces total disk space usage and is faster.

See the [documentation for `nf-core download`](/docs/nf-core-tools/pipelines/download) for more information.

### Reference genomes

Some pipelines require reference genomes and have built-in integration with AWS-iGenomes.
If you wish to use these references, you must also download and transfer them to your offline system.

Follow the [reference genomes documentation](/docs/usage/reference_genomes/reference_genomes.md) to configure the base path for the references.

### Bytesize talk

Here is a recent bytesize talk explaining the necessary steps to run pipelines offline.

<!-- markdownlint-disable -->
<iframe width="560" height="315" src="https://www.youtube.com/embed/N1rRr4J0Lps" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
<!-- markdownlint-restore -->
