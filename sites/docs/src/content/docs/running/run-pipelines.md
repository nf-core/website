---
title: Running pipelines
subtitle: Run nf-core pipelines
shortTitle: Running pipelines
---

With Nextflow and your software dependency manager installed, you are ready to run nf-core pipelines.
This guide covers the essential commands and patterns for running any nf-core pipeline.

:::tip
For a hands-on introduction, see the [Run your first pipeline](../get_started/run-your-first-pipeline.md) tutorial using [nf-core/demo](https://github.com/nf-core/demo).
:::

## Basic command structure

All nf-core pipelines follow a consistent command structure:

```bash
nextflow run nf-core/<pipeline> -r <version> -profile <config-name> <pipeline-parameters...>
```

- `nf-core/<pipeline>`: The pipeline name (for example `nf-core/rnaseq` or `nf-core/sarek`)
- `-r <version`>: the version, tag, or commit of the pipeline to execute (for example `1.2.0`, `b52fa53`). Optional but recommended
- `-profile <config-name>`: Configuration profile(s) for software dependencies and execution environment
- `<pipeline-parameters...>`: Pipeline-specific parameters and Nextflow options

## Common operations

The following sections cover common operations for running nf-core pipelines.

### Testing with test data

All nf-core pipelines include a `test` profile with small datasets for verifying your setup:

```bash
nextflow run nf-core/<pipeline> -profile test,<config> --outdir results
```

- Replace `<pipeline>` with the pipeline name
- Replace `<config>` with your software dependency manager

:::tip
Common software dependency and compute environment profiles include:

- `docker`: Uses Docker containers
- `singularity`: Uses Singularity containers
- `conda`: Uses Conda environments

Combine profiles with commas: `-profile test,docker`.
For a full list of supported options, see [Software dependencies](./environment_setup/software-dependencies.md).

<!-- TODO software-dependencies is a dead link -->

:::

:::warning
nf-core test profiles require an internet connection.
:::

Results are written to the `results/` directory, organized by analysis step. All pipelines include a `pipeline_info/` subdirectory with execution logs and metadata.

### Running with your own data

Most nf-core pipelines use a samplesheet to specify input files:

1. Create a CSV file with your input data. The required columns vary by pipeline—check the pipeline documentation for the specific format.

   Example samplesheet structure:

   ```csv
   sample,fastq_1,fastq_2
   SAMPLE1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
   SAMPLE2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
   ```

2. Run the pipeline with your samplesheet:

   ```bash
   nextflow run nf-core/<pipeline> \
     -profile <config> \
     --input samplesheet.csv \
     --outdir results
   ```

:::tip
Each pipeline's documentation includes:

- Required samplesheet format and columns
- Example samplesheets
- Details about required and optional parameters

Consult the pipeline's usage and parameter documentation before preparing your data.
:::

### Using parameter files

Parameter files let you specify pipeline options in a structured format rather than on the command line:

1. Create a JSON or YAML file with your parameters:

   ```json
   {
     "input": "samplesheet.csv",
     "outdir": "results",
     "max_memory": "128.GB",
     "max_cpus": 16
   }
   ```

2. Run the pipeline with your parameter file:

   ```bash
   nextflow run nf-core/<pipeline> \
     -profile <config> \
     -params-file params.json
   ```

### Configuring pipelines

For complex execution environments or custom compute resource requirements (for example, memory or CPU adjustments), use configuration files to override default settings.

#### Using configuration files

Create a custom configuration file to specify execution settings:

1. Create a file (for example, `custom.config`) with your configuration:

   ```groovy
   // Resource requests for all processes
   process {
     cpus   = 8
     memory = 32.GB
     time   = 24.h
   }

   // Executor settings
   executor {
     name = 'slurm'
     queueSize = 50
   }
   ```

2. Apply your configuration with the `-c` flag:

   ```bash
   nextflow run nf-core/<pipeline> \
     -profile <config> \
     -c custom.config \
     --input samplesheet.csv \
     --outdir results
   ```

#### Common configuration scenarios

Common configuration scenarios include:

- Adjust resource requests

  ```groovy
  process {
    withName: 'PROCESS_NAME' {
      cpus   = 16
      memory = 64.GB
    }
  }
  ```

- Limit maximum resource requests

  ```groovy
  process {
    resourceLimits = [
        memory: 1000.GB,
        cpus: 144,
        time: 30.d,
    ]
  }
  ```

- Configure cluster execution

  ```groovy
  process {
    executor = 'slurm'
    queue    = 'high-priority'
    clusterOptions = '--account=myproject'
  }
  ```

- Set custom paths

  ```groovy
  params {
    custom_config_base = '/path/to/configs'
  }
  ```

:::tip
Configuration files are applied in order and later files override earlier settings. Multiple `-c` flags can be used to layer configurations.
:::

For detailed configuration options, see the [Configuration guide](../running/configuration/overview.md).

<!-- TODO Fix the configuration guide link -->

### Resuming runs

Add `-resume` to continue from where a pipeline stopped (due to errors or cancellation):

```bash
nextflow run nf-core/<pipeline> \
  -profile <config> \
  --input samplesheet.csv \
  --outdir results \
  -resume
```

Nextflow caches completed steps and only re-runs processes with changed inputs, saving time and resources.

:::tip
Always use `-resume` when re-running a pipeline after fixing errors or updating parameters. Nextflow automatically determines which steps need to re-execute.
:::

## Next steps

- [Configure pipelines](../running/configuration/overview.md) to adjust resource requirements and execution settings <!-- TODO deadlink -->
- Browse the [pipeline catalog](https://nf-co.re/pipelines) to find workflows for your research
- Join the [nf-core Slack](https://nf-co.re/join/slack) community for support

<!-- TODO: Add links, if any extra seem applicable -->
