---
title: Configuration options
subtitle: Configure pipelines
shortTitle: Configuration options
---

You can configure pipelines using three approaches:

- [Default configuration profiles](#default-configuration-profiles)
- [Shared nf-core/configs](#shared-nf-coreconfigs)
- [Custom configuration files](#custom-configuration-files)

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

- You need pipeline-specific resource limits
- Running on unique infrastructure
- You are the only user of the pipeline

### Default configuration profiles

Each nf-core pipeline includes default resource requirements in `config/base.config`.
The pipeline loads this base configuration first, then overwrites it with any subsequent configuration layers.

Enable configuration profiles using the `-profile` command line flag.
You can specify multiple profiles in a comma-separated list (e.g., `-profile test,docker`).

:::note
Order matters.
Profiles load in sequence.
Later profiles overwrite earlier ones.
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

Without a specified profile, the pipeline runs locally and expects all software to be installed and available on the `$PATH`.
This approach is not recommended.

Each pipeline includes `test` and `test_full` profiles.
These run the workflow with minimal or full-size public datasets for automated CI tests.
You can use them to test the pipeline on your infrastructure before running your own data.

### Shared nf-core/configs

If you work on a shared system (for example, an HPC cluster or institutional server), use a configuration profile from [nf-core/configs](https://github.com/nf-core/configs).
All nf-core pipelines load these shared profiles at run time.

Check if your system has a profile at [https://github.com/nf-core/configs](https://github.com/nf-core/configs).
If not, follow the repository instructions or the tutorial to add your cluster.

See [Run configuration](https://training.nextflow.io/latest/nextflow_run/03_config/) for a guided configuration tutorial.

:::note
All nf-core pipelines have a `conf/base.config` where default computational resource requests are specified.
:::

### Custom configuration files

You can add additional configuration files to customize your pipeline execution.

When you launch a pipeline script, Nextflow detects configuration files from multiple sources and applies them in the following order (from lowest to highest priority):

1. `$NXF_HOME/config` (defaults to `$HOME/.nextflow/config`)
1. `nextflow.config` in the project directory
1. `nextflow.config` in the launch directory
1. Config files specified with `-c <config-files>`

Nextflow loads configuration **parameters** sequentially and overwrites previous values.
The loading order is:

1. Pipeline defaults
1. User's home directory
1. Working directory
1. Each `-c` file in the order you specify
1. Command line parameters (`--<parameter>`)

:::warning
Parameters set in your `custom.config` files will not override defaults in `nextflow.config`.
Use `-params-file` with YAML or JSON format instead.
:::

:::tip
Generate a parameters file using the **Launch** button on the [nf-co.re website](https://nf-co.re/launch).
:::

## Additional resources

For more information about configuration syntax and parameters, see:

- [Nextflow config](https://www.nextflow.io/docs/latest/config.html)
- [Institutional profiles](../../developing/institutional-profiles/overview.md)
