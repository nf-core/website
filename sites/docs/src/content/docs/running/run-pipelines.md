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
nextflow run nf-core/<pipeline> -r <version> -profile <config-name> <pipeline-parameters>
```

- `nf-core/<pipeline>`: The pipeline name (for example, `nf-core/rnaseq` or `nf-core/sarek`).
- `-r <version>`: The version, tag, or commit of the pipeline to execute (for example, `1.2.0`, `b52fa53`). Optional but recommended.
- `-profile <config-name>`: Configuration profile(s) for software dependencies and execution environment.
- `<pipeline-parameters>`: Pipeline-specific parameters and Nextflow options.

:::note
Nextflow options such as `-r`, `-profile`, and `-resume` use a single hyphen.
Pipeline parameters such as `--input` and `--outdir` use a double hyphen.
:::

## Common operations

The following sections cover common operations for running nf-core pipelines.

### Run with test data

All nf-core pipelines include a `test` profile with small datasets for verifying your setup:

```bash
nextflow run nf-core/<pipeline> -profile test,<config> --outdir results
```

- Replace `<pipeline>` with the pipeline name.
- Replace `<config>` with your software dependency manager.

:::tip
Common software dependency and compute environment profiles include:

- `docker`: Uses Docker containers
- `singularity`: Uses Singularity containers
- `conda`: Uses Conda environments

Combine profiles with commas: `-profile test,docker`.
For a full list of supported options, see [Software dependencies](../get_started/environment_setup/software-dependencies).

:::

:::warning
nf-core test profiles require an internet connection.
:::

Results are written to the `results/` directory, organized by analysis step. All pipelines include a `pipeline_info/` subdirectory with execution logs and metadata.

### Run with your own data

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

### Specify a pipeline version

Specify the pipeline version with `-r` when you run a pipeline on your data.
This ensures reproducibility of your runs, for example when publishing your analysis.
:::warning
The first time you pin a version it keeps the pipeline code and software versions the same for every run, even after the pipeline changes on GitHub.

Make sure to explicitly update the version when you want to use a new version of the pipeline.
:::

1. Find the version you want on the pipeline's nf-core page under the 'Releases' tab at https://nf-co.re/<pipeline>/releases_stats/, or on the pipeline's GitHub releases page, at `https://github.com/nf-core/<pipeline>/releases`. Use the numeric tag only, for example `1.3.1`.

2. Add the version to your command with `-r`:

   ```bash
   nextflow run nf-core/<pipeline> \
     -r 1.3.1 \
     -profile <config> \
     --input samplesheet.csv \
     --outdir results
   ```

To switch to another version, change the value after `-r`.

The pipeline records the version in the execution reports in the `pipeline_info/` subdirectory of your results, and at the bottom of the MultiQC report for pipelines that use MultiQC.
Use these to check which version produced a set of results.

### Use parameter files

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

:::tip
Parameter files make runs easier to reproduce and share, for example as supplementary material for a publication.
Before you share a parameter file, remove any paths that are specific to your cluster and any institutional profiles.
:::

### Configure pipelines

For complex execution environments or custom compute resource requirements (for example, memory or CPU adjustments), use configuration files to override default settings.

#### Use configuration files

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

For detailed configuration options, see the [Configuration guide](./configuration/configuration-options).

### Resume runs

Add `-resume` to continue from where a pipeline stopped (due to errors or cancellation):

```bash
nextflow run nf-core/<pipeline> \
  -profile <config> \
  --input samplesheet.csv \
  --outdir results \
  -resume
```

Nextflow caches completed steps and only re-runs processes with changed inputs, saving time and resources.
Nextflow treats an input as unchanged only if both the file name and the file contents are identical.
For more information, see the [Demystifying Nextflow resume](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html) blog post.

:::tip
Always use `-resume` when re-running a pipeline after fixing errors or updating parameters. Nextflow automatically determines which steps need to re-execute.
:::

To resume an earlier run instead of the most recent one, add its run name:

```bash
nextflow run nf-core/<pipeline> \
  -profile <config> \
  --input samplesheet.csv \
  --outdir results \
  -resume <run-name>
```

Run `nextflow log` to list the names of previous runs.

### Run in the background

Nextflow submits and supervises pipeline jobs, and must keep running until the pipeline finishes.
If you exit a terminal while Nextflow runs in the foreground, the pipeline stops.

To keep a pipeline running after you exit a terminal, use one of these options:

- Add `-bg` to launch Nextflow in the background, detached from your terminal. Nextflow writes the log output to a file instead of the terminal.

  ```bash
  nextflow run nf-core/<pipeline> \
    -profile <config> \
    --input samplesheet.csv \
    --outdir results \
    -bg
  ```

- Use `screen`, `tmux`, or a similar tool to create a detached session that you can reattach to later.
- On some high-performance computing (HPC) systems, it is required to submit Nextflow itself as a job to the cluster scheduler. Nextflow then submits the pipeline jobs from inside that job.

### Update pipelines

The first time you run a pipeline, Nextflow pulls the code from GitHub and caches a copy in `~/.nextflow/assets/`.
Later runs use the cached copy, even if the pipeline has changed on GitHub since.

To run the latest version of a pipeline, update the cached copy regularly:

```bash
nextflow pull nf-core/<pipeline>
```

:::tip
If an update does not take effect, see [Update not applying](./troubleshooting#update-not-applying).
:::

### Notifications

- Specify `--email your@email.com` to receive emails when your pipeline completes (requires Nextflow [mail and notification configuration](https://docs.seqera.io/nextflow/notifications#mail)).
- To receive notifications from your pipeline (e.g., on completion or failure), see the [nf-slack plugin](https://github.com/seqeralabs/nf-slack) for Slack or the [nf-teams plugin](https://github.com/nvnieuwk/nf-teams) for Microsoft Teams.

## Next steps

For more resources to help you run pipelines:

- Learn how to [Configure pipelines](./configuration/configuration-options)
- Browse the [pipeline catalog](https://nf-co.re/pipelines) to find workflows for your research
- Join the [nf-core Slack](https://nf-co.re/join/slack) for community support
