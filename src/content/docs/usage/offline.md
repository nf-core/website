---
title: Running offline
subtitle: How to use nf-core pipelines without an internet connection.
---

# Introduction

Nextflow supports fetching nearly everything it needs to run a pipeline over the web automatically: pipeline code, software requirements, reference genomes, and even remote data sources.

Nextflow can run your analysis on a system that has no internet connection.
There are just a few extra steps required to get everything you need available locally.

You will need to fetch everything you need on a system that _does_ have an internet connection (typically your personal computer).
Then, transfer these to your offline system using the methods you have available.

You will need three things: a working version of [Nextflow](#nextflow), the [pipeline assets](#pipeline-code), and any required [reference genomes](#reference-genomes).

## Nextflow

First of all, you need to have Nextflow installed on your system.
We do this by installing it locally on a machine that _does_ have an internet connection, and then transferring to the offline system.

- Start by [installing Nextflow locally](https://nextflow.io/docs/latest/getstarted.html#installation).
  :::warning
  do _not_ use the `-all` package, as this does not allow the use of custom plugins.
  :::
- Kick off a pipeline locally so that Nextflow fetches the required plugins. It does not need to run to completion.
- Copy the Nextflow binary and `$HOME/.nextflow` folder to your offline environment.
- In your Nextflow configuration file, specify each plugin that you downloaded, both name and version, including default plugins. This will prevent Nextflow from trying to download newer versions of plugins.
- Add the following environment variable in your `~/.bashrc` file:
  ```bash title=".bashrc"
  export NXF_OFFLINE='true'
  ```

## Pipeline code

To run a pipeline offline, you need the pipeline code, the software dependencies, and the shared nf-core/configs configuration profiles.
We have created a helper tool as part of the _nf-core_ package to automate this for you.

On a computer with an internet connection, run `nf-core download <pipeline>` to download the pipeline and config profiles.
Add the argument `--container singularity` to also fetch the singularity container(s).

The pipeline and requirements will be downloaded, configured with their relative paths, and packaged into a `.tar.gz` file by default.
This can then be transferred to your offline system and unpacked.

Inside, you will see directories called `workflow` (the pipeline files), `config` (a copy of [nf-core/configs](https://github.com/nf-core/configs)), and (if you used `--container singularity`) a directory called `singularity`.
The pipeline code is adjusted by the download tool to expect these relative paths, so as long as you keep them together it should work as is.

To run the pipeline, simply use `nextflow run <download_directory>/workflow [pipeline flags]`.

### Shared storage

If you are downloading _directly_ to the offline storage (eg. a head node with internet access whilst compute nodes are offline), you can use the `--singularity-cache-only` option for `nf-core download` and set the `$NXF_SINGULARITY_CACHEDIR` environment variable.
This downloads the singularity images to the `$NXF_SINGULARITY_CACHEDIR` folder and does not copy them into the target downloaded pipeline folder.
This reduces total disk space usage and is faster.

For more information, see the [documentation for `nf-core download`](https://nf-co.re/tools#downloading-pipelines-for-offline-use).

## Reference genomes

Some pipelines require reference genomes and have built-in integration with AWS-iGenomes.
If you wish to use these references, you must download and transfer them to your offline cluster.
Once transferred, follow the [reference genomes documentation](reference_genomes.md) to configure the base path for the references.

## Bytesize talk

Here is a bytesize talk explaining the necessary steps to run pipelines offline.

<!-- markdownlint-disable -->
<iframe width="560" height="315" src="https://www.youtube.com/embed/N1rRr4J0Lps" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
<!-- markdownlint-restore -->
