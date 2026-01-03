---
title: Configuration
subtitle: Learn about nf-core pipeline configuration
shortTitle: Overview
weight: 1
---

Nextflow runs on HPC execution schedulers ([Slurm](https://slurm.schedmd.com/quickstart.html), [SGE](https://docs.oracle.com/cd/E19680-01/html/821-1541/ciagcgha.html#scrolltoc), and [others](https://www.nextflow.io/docs/latest/executor.html)) and cloud platforms ([AWS Batch](https://aws.amazon.com/batch/), [Google Cloud](https://cloud.google.com/), and [others](https://www.nextflow.io/docs/latest/executor.html)).

Nextflow also supports container engines ([Docker](https://www.docker.com/), [Singularity](https://sylabs.io/), and [others](https://www.nextflow.io/docs/latest/container.html)) and dependency managers ([Conda](https://docs.conda.io/en/latest/) and [Spack](https://spack.readthedocs.io/en/latest/)) for software deployment.

To run nf-core pipelines on your system, install your dependency management software (see [Installation](/docs/usage/getting_started/installation)) and configure Nextflow.

:::note{title="Configuration vs parameters"}
Configuration controls how Nextflow runs (executor, resources, containers). Parameters control what the pipeline does with your data (`--input`, `--outdir`).

This section covers configuration. For pipeline-specific parameters, see the pipeline documentation.
:::

## Configuration options

You can configure pipelines using three approaches:

1. [Default pipeline configuration profiles](#default-configuration-profiles)
2. [Shared nf-core/configs configuration profiles](#shared-nf-coreconfigs)
3. [Custom configuration files](#custom-configuration-files)

:::warning{title="Do not edit the pipeline code to configure nf-core pipelines"}
Editing pipeline defaults prevents you from updating to newer versions without overwriting your changes.
This breaks reproducibility and moves your workflow away from the canonical pipeline.
:::

## Choosing your configuration approach

Use default profiles (`-profile docker,test`) when:

- Running quick tests
- Using standard container engines
- Running on your local machine

Use shared nf-core/configs when:

- Working on a shared HPC cluster
- Your institution already has a profile
- Multiple users run pipelines on the same system

Use custom configuration files when:

- You need specific resource limits
- Running on unique infrastructure
- You are the only user of the pipeline

### Default configuration profiles

Each nf-core pipeline includes default resource requirements in `config/base.config`.
The pipeline loads this base configuration first, then overwrites it with any subsequent configuration layers.

Enable configuration profiles using the `-profile` command line flag.
You can specify multiple profiles in a comma-separated list (e.g., `-profile test,docker`).

:::note
Order matters. Profiles load in sequence. Later profiles overwrite earlier ones.
:::
nf-core provides these basic profiles for container engines:

- `docker`: Uses [Docker](http://docker.com/) and pulls software from quay.io
- `apptainer`: Uses [Singularity](https://apptainer.org/) and pulls software from quay.io
- `podman`: Uses [Podman](https://podman.io/)
- `shifter`: Uses [Shifter](https://www.nersc.gov/research-and-development/user-defined-images/)
- `charliecloud`: Uses [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`: Uses [Conda](https://conda.io/docs/) and pulls most software from [Bioconda](https://bioconda.github.io/)

:::note
Use Conda only as a last resort (that is, when you cannot run the pipeline with Docker or Singularity).
:::

Without a specified profile, the pipeline runs locally and expects all software to be installed and available on the `PATH`.
This approach is not recommended.

Each pipeline includes `test` and `test_full` profiles.
These run the workflow with minimal or full-size public datasets for automated CI tests.
You can use them to test the pipeline on your infrastructure before running your own data.

### Shared nf-core/configs

If you work on a shared system (for example, an HPC cluster or institutional server), use a configuration profile from [nf-core/configs](https://github.com/nf-core/configs).
All nf-core pipelines load these shared profiles at run time.

Check if your system has a profile at [https://github.com/nf-core/configs](https://github.com/nf-core/configs).
If not, follow the repository instructions or the tutorial to add your cluster.

<!-- TODO: Add link to tutorial -->

### Custom configuration files

If you run the pipeline alone, create a local configuration file.
Nextflow searches for configuration files in three locations:

1. User's home directory: `~/.nextflow/config`
2. Analysis working directory: `nextflow.config`
3. Custom path on the command line: `-c path/to/config` (you can specify multiple files)

Nextflow loads configuration parameters sequentially and overwrites previous values.
The loading order is:

1. Pipeline defaults
2. User's home directory
3. Working directory
4. Each `-c` file in the order you specify
5. Command line parameters (`--<parameter>`)

:::warning
Parameters in `custom.config` files will not override defaults in `nextflow.config`.
Use `-params-file` with YAML or JSON format instead.
:::

:::tip
Generate a parameters file using the **Launch** button on the [nf-co.re website](https://nf-co.re/launch).
:::

## Additional resources

For more information about configuration syntax and parameters, see:

- [Nextflow config](https://www.nextflow.io/docs/latest/config.html)

<!-- TODO: Add links, if any -->
