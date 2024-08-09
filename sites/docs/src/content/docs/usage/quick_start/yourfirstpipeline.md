---
title: Your first pipeline
subtitle: Rune your first Nextflow and nf-core pipelines.
shortTitle: Your first pipeline
weight: 3
---

## Your first pipeline

Nextflow works best with an active internet connection, as it is able to fetch all pipeline requirements. See [Running offline](/docs/usage/configuration/runningoffline) if you need to run Nextflow pipelines without an internet connection.

### Hello world!

Run Nextflow’s [Hello world!](https://github.com/nextflow-io/hello) pipeline example to confirm Nextflow is installed and that you have an active internet connection:

```bash
nextflow run hello
```

If everything has been installed properly you should see the pipeline launch and run to completion.

### nf-core pipelines

Each nf-core pipeline is different and will have different run requirements.

Before running a pipeline with your own data it is best to test it is working using a small dataset that is known to work.

nf-core pipelines come packed with a `test` profile that uses small data files from the [nf-core/test-datasets](https://github.com/nf-core/test-datasets) repository. The `test` profile is extensively tested and should run without extensive configuration. nf-core pipelines also come packed with directives for containers and environments that can be flexibly enabled using software dependency profiles (e.g., `docker`, `singularity`, and `conda`). By using a profile for software dependencies an nf-core pipeline should run with all required tools.

Run the `nf-core/rnaseq` pipeline using the `test` profile and a software profile (e.g., `docker`) to confirm Nextflow is installed and that you have an active internet connection:

```bash
nextflow run nf-core/rnaseq -profile test,docker --outdir results
```

If everything has been installed properly you should see the pipeline launch and run to completion.

You are now ready to start adding your own data and configure the pipeline to better suit your requirements.

See the [Configuration](/docs/usage/configuration/introduction) section for more information about how you can customize the execution of a pipeline.
