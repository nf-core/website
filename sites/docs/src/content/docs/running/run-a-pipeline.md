---
title: Running pipelines
subtitle: Run an nf-core pipeline
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
nextflow run nf-core/<PIPELINE> -profile <CONFIG> [OPTIONS]
```

- `nf-core/<PIPELINE>`: The pipeline name (e.g., `nf-core/rnaseq`, `nf-core/sarek`)
- `-profile <CONFIG>`: Configuration profile(s) for software dependencies and execution environment
- `[OPTIONS]`: Pipeline-specific parameters and Nextflow options

## Common operations

The following sections cover common operations for running nf-core pipelines.

### Testing with test data

All nf-core pipelines include a `test` profile with small datasets for verifying your setup:

```bash
nextflow run nf-core/<PIPELINE> -profile test,<CONFIG> --outdir results
```

- Replace `<PIPELINE>` with the pipeline name
- Replace `<CONFIG>` with your software dependency manager

:::tip
Common software dependency and compute environment profiles include:

- `docker`: Uses Docker containers
- `singularity`: Uses Singularity containers
- `conda`: Uses Conda environments

Combine profiles with commas: `-profile test,docker`.
For a full list of supported options, see [Software dependencies](./environment_setup/software-dependencies.md).
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
   nextflow run nf-core/<PIPELINE> \
     -profile <CONFIG> \
     --input samplesheet.csv \
     --outdir results
   ```

:::tip
Each pipeline's documentation includes:

- Required samplesheet format and columns
- Example samplesheets
- Details about required and optional parameters

Consult the pipeline's usage documentation before preparing your data.
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
   nextflow run nf-core/<PIPELINE> \
     -profile <CONFIG> \
     -params-file params.json
   ```

### Configuring pipelines

For complex execution environments or custom resource requirements, use configuration files to override default settings.

#### Using configuration files

Create a custom configuration file to specify execution settings:

1. Create a file (e.g., `custom.config`) with your configuration:

   ```groovy
   // Resource limits
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
   nextflow run nf-core/<PIPELINE> \
     -profile <CONFIG> \
     -c custom.config \
     --input samplesheet.csv \
     --outdir results
   ```

#### Common configuration scenarios

Common configuration scenarios include:

- **Adjust resource limits**

  ```groovy
  process {
    withName: 'PROCESS_NAME' {
      cpus   = 16
      memory = 64.GB
    }
  }
  ```

- **Configure cluster execution**

  ```groovy
  process {
    executor = 'slurm'
    queue    = 'high-priority'
    clusterOptions = '--account=myproject'
  }
  ```

- **Set custom paths**

  ```groovy
  params {
    custom_config_base = '/path/to/configs'
    max_memory = '128.GB'
    max_cpus = 32
  }
  ```

:::tip
Configuration files are applied in order and later files override earlier settings. Multiple `-c` flags can be used to layer configurations.
:::

For detailed configuration options, see the [Configuration guide](../running/configuration/overview.md).

### Resuming runs

Add `-resume` to continue from where a pipeline stopped (due to errors or cancellation):

```bash
nextflow run nf-core/<PIPELINE> \
  -profile <CONFIG> \
  --input samplesheet.csv \
  --outdir results \
  -resume
```

Nextflow caches completed steps and only re-runs processes with changed inputs, saving time and resources.

:::tip
Always use `-resume` when re-running a pipeline after fixing errors or updating parameters. Nextflow automatically determines which steps need to re-execute.
:::

### Specifying pipeline versions

Use `-r` to specify a pipeline version, tag, or commit:

```bash
nextflow run nf-core/<PIPELINE> \
  -r 1.0.2 \
  -profile <CONFIG> \
  --input samplesheet.csv \
  --outdir results
```

## Next steps

- [Configure pipelines](../running/configuration/overview.md) to adjust resource requirements and execution settings
- Browse the [pipeline catalog](https://nf-co.re/pipelines) to find workflows for your research
- Join the [nf-core Slack](https://nf-co.re/join/slack) community for support
