---
title: Running offline
subtitle: Using nf-core pipelines without an internet connection.
---

Nextflow supports fetching nearly everything it needs to run a pipeline over the web automatically: pipeline code, software requirements, reference genomes and even remote data sources.

If you need to run your analysis on a system that has no internet connection, don't panic!
There are just a few extra steps required to get everything you need available locally.

## Nextflow

First of all, you need to have Nextflow installed on your system.
Go to the Nextflow releases page on GitHub: [https://github.com/nextflow-io/nextflow/releases](https://github.com/nextflow-io/nextflow/releases).
Each release has a dropdown with associated _Assets_.
One of these should have the suffix `-all`, _e.g._ `nextflow-19.10.0-all`.
Download this file and transfer to your offline system.
Run it to install Nextflow (it is a very large _bash_ file).

Once installed, you can stop nextflow from looking for updates online by adding the following environment variable in your `~/.bashrc` file:

```bash
export NXF_OFFLINE='TRUE'
```

## Pipeline code

To run a pipeline offline you need the pipeline code, the software requirements and the shared nf-core/configs configuration profiles.
To help with this process, we have created a helper tool as part of the _nf-core_ package to automate this for you.

On a computer with an internet connection, run `nf-core download <pipeline>` to download the pipeline and config profiles.
Add the argument `--singularity` to also fetch the singularity container (requires Singularity to be installed on your system).

The pipeline and requirements will be downloaded, configured with their relative paths and packaged in to a `.tar.gz` file by default. This can then be transferred to your offline system and unpacked.

To run the pipeline, simply do `nextflow run <directory>`

For more information about the nf-core helper tools, see the [documentation](https://nf-co.re/tools#downloading-pipelines-for-offline-use).

## Reference genomes

Some pipelines require reference genomes and have built-in integration of AWS-iGenomes.
If you wish to use these references you must download them and transfer to your offline cluster.
Once transferred, follow the [reference genomes documentation](reference_genomes.md) to configure the base path for the references.
