---
title: Configuring pipelines for your system
subtitle: Learn how to configure pipelines for your system
shortTitle: Configuring pipelines for your system
weight: 2
---

nf-core pipelines work on different computational infrastructures, from laptops to HPC clusters and cloud platforms.
This page shows you how to configure pipelines to match your system's capabilities.

## Understand workflow resources

The base configuration of nf-core pipelines defines default resource allocations for each workflow step (e.g., in the [`base.config`](https://github.com/nf-core/rnaseq/blob/master/conf/base.config) file).

These default values are generous to accommodate diverse workloads across different users.
Your jobs might receive more resources than needed, which can reduce system efficiency.
You might also want to increase resources for specific tasks to maximise speed.
Consider increasing resources if a pipeline step fails with a `Command exit status` of `137`.

Pipelines configure tools to use available resources when possible (e.g., with `-p ${task.cpus}`), where `${task.cpus}` is dynamically set from the pipeline configuration.
Not all tools support dynamic resource configuration.

Most process resources use process labels, as shown in this base configuration example:

```groovy
process {
  resourceLimits = [
    cpus: 32,
    memory: 256.GB,
    time: 24.h
  ]
  withLabel:process_low {
    cpus   = { 2 * task.attempt }
    memory = { 14.GB * task.attempt }
    time   = { 6.h  * task.attempt }
  }
  withLabel:process_medium {
    cpus   = { 6  * task.attempt }
    memory = { 42.GB * task.attempt }
    time   = { 8.h * task.attempt }
  }
  withLabel:process_high {
    cpus   = { 12 * task.attempt }
    memory = { 84.GB * task.attempt }
    time   = { 10.h * task.attempt }
  }
}
```

The `resourceLimits` list sets the absolute maximum resources any pipeline job can request (typically matching your machine's maximum available resources).
The label blocks define the initial default resources each pipeline job requests.
When a job runs out of memory, most nf-core pipelines retry the job and increase the resource request up to the `resourceLimits` maximum.

### Customize process resources

:::tip
Copy only the labels you want to change into your custom configuration file, not all labels.
:::

To set a fixed memory allocation for all large tasks across most nf-core pipelines (without increases during retries), add this to your custom configuration file:

```groovy
process {
  withLabel:process_high {
    memory = 200.GB
  }
}
```

You can target a specific process (job) name instead of a label using `withName`.
Find process names in your console log when the pipeline runs.
For example:

```groovy
process {
  withName: STAR_ALIGN {
    cpus = { 32 * task.attempt }
  }
}
```

When a pipeline uses a tool multiple times in the workflow, specify the complete execution path of the module:

```groovy
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

:::info
If you receive a warning about an unrecognised process selector, check that you specified the process name correctly.
:::

For more information, see the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#process-selectors).

After writing your [configuration file](#custom-configuration-files), supply it to your pipeline command with `-c`.

:::warning
Check your syntax carefully. Setting memory to `200` allocates 200 bytes of memory.

Use quotes with a space or no quotes with a dot: `"200 GB"` or `200.GB`.

See the Nextflow documentation for [memory](https://www.nextflow.io/docs/latest/process.html#memory), [cpus](https://www.nextflow.io/docs/latest/process.html#cpus), and [time](https://www.nextflow.io/docs/latest/process.html#time).
:::

If the pipeline defaults need adjustment, contact the pipeline developers on Slack in the pipeline channel or submit a GitHub issue on the pipeline repository.

## Change your executor

Nextflow pipelines run in local mode by default, executing jobs on the same system where Nextflow runs.
Most users need to specify an executor to tell Nextflow how to submit jobs to a job scheduler (e.g., SGE, LSF, Slurm, PBS, or AWS Batch).

You can configure the executor in shared configuration profiles or in custom configuration files.
For more information about executors, see [Executors](https://www.nextflow.io/docs/latest/executor.html#executor-page).

## Set max resources

Pipeline runs can fail when a step requests more resources than your system has available.

To prevent these failures, configure Nextflow to cap resource requests using the `resourceLimits` list in your configuration file.
These limits should match the maximum resources of your machine or node.

For example, add this to your Nextflow configuration:

```groovy
process {
  resourceLimits = [
    cpus: 32,
    memory: 256.GB,
    time: 24.h
  ]
}
```

When a job exceeds the default memory request, Nextflow retries the job with increased memory.
The memory increases with each retry until the job completes or reaches the `256.GB` limit.

These parameters cap resource requests to prevent Nextflow from submitting jobs that exceed your system's capabilities.

Specifying resource limits does not increase the resources available to pipeline tasks.
See [Tuning workflow resources](#tuning-workflow-resources) for more information.

:::note{collapse title="Note on older nf-core pipelines"}

For pipelines generated with the nf-core template before v3.0.0, or when running with Nextflow versions earlier than 24.04.0, you may need to use a different syntax to prevent resources from exceeding a maximum limit.

Pipeline runs can fail when a step requests more resources than your system has available.

To avoid these failures, all nf-core pipelines [check](https://github.com/nf-core/tools/blob/99961bedab1518f592668727a4d692c4ddf3c336/nf_core/pipeline-template/nextflow.config#L206-L237) pipeline-step resource requests against parameters called `--max_cpus`, `--max_memory`, and `--max_time`.
These should represent the maximum possible resources of a machine or node.

These parameters act as a cap to prevent Nextflow from submitting a single job requesting more resources than your system has available.

Increasing these values from the defaults does not increase the resources available to the pipeline tasks.
See [Tuning workflow resources](#tuning-workflow-resources) for more information.

Most pipelines automatically restart jobs that fail due to lack of resources with doubled requests.
These caps keep those requests from exceeding system limits and crashing the entire pipeline run.
When a job exceeds the process-specific default resources and is retried, only resource requests (CPU, memory, or time) that have not yet reached the value set with `--max_<resource>` will be increased during the retry.

The `--max_<resource>` parameters represent the maximum for a single pipeline job, not the total sum of resource usage at any given time.

:::

## Customize Docker registries

Most pipelines use `quay.io` as the default Docker registry for Docker and Podman images.
When you specify a Docker container without a full URI, Nextflow pulls the image from `quay.io`.

For example, this container specification:

- `biocontainers/fastqc:0.11.7--4`

Pulls from `quay.io`, resulting in the full URI:

- `quay.io/biocontainers/fastqc:0.11.7--4`

If you specify a different `docker.registry` value, Nextflow uses that registry instead.
For example, if you set `docker.registry = 'myregistry.com'`, the image pulls from:

- `myregistry.com/biocontainers/fastqc:0.11.7--4`

When you specify a full URI in the container specification, Nextflow ignores the `docker.registry` setting and pulls exactly as specified:

- `docker.io/biocontainers/fastqc:v0.11.9_cv8`

## Update tool versions

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of nf-core pipelines uses one container or Conda environment per process, which simplifies software dependency maintenance and updates.

To use a different version of a tool, identify the `process` name and override the Nextflow `container` or `conda` definition for that process using the `withName` declaration.

For example, the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline uses Pangolin, which updates its internal COVID-19 lineage database frequently.
Re-releasing nf-core/viralrecon with every Pangolin update would be impractical.

You can override the default container by creating a custom configuration file and passing it with `-c custom.config`.

1. Check the default version in the module file in the `modules/nf-core/` directory of the pipeline (e.g., see [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19))
1. Find the latest Biocontainer version on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags) for Docker or [Galaxy Project](https://depot.galaxyproject.org/singularity/) for Singularity
   - The container version tag is identical for both systems but must include the build ID (e.g., `--pyhdfd78af_1`)
1. Create the custom configuration:

   - For Docker:

     ```groovy
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.1.17--pyhdfd78af_1'
         }
     }
     ```

   - For Singularity:

     ```groovy
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.1.17--pyhdfd78af_1'
         }
     }
     ```

   - For Conda (check against [bioconda](https://bioconda.github.io), which does not include the build tag):

     ```groovy
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.1.17'
         }
     }
     ```

:::warning
Pipeline developers provide no warranty when you update containers.
Major changes in the container tool may break the pipeline.
:::

:::warning
Tool developers sometimes change version reporting between updates.
Container updates may break version reporting within the pipeline and create missing values in MultiQC version tables.
:::

## Modifying tool arguments

You might need to understand which arguments a pipeline uses for a specific tool, or modify these arguments for your analysis.

Check the tool parameter descriptions in the pipeline documentation.
For example, see the **Help** button next to [this parameter](https://nf-co.re/funcscan/2.0.0/parameters#annotation_bakta_mincontiglen) in [nf-core/funcscan](https://nf-co.re/funcscan).

### Find tool arguments

If the parameter description doesn't list the arguments, check these two locations:

Most arguments (mandatory and optional) are defined in `conf/modules.conf` under the `ext.args` entry.
See the default arguments for the `SORTMERNA` step in nf-core/rnaseq [here](https://github.com/nf-core/rnaseq/blob/6e1e448f535ccf34d11cc691bb241cfd6e60a647/conf/modules.config#LL299).

:::info
The pipeline inserts `ext.args` arguments into the module using the `$args` variable in the module's Bash code.
:::

Some modules have mandatory arguments specified directly in the module code.
Find these in the `script` section of the module file at `modules/<nf-core/local>/<tool>/main.nf` in the pipeline's GitHub repository.

### Customise tool arguments

To modify tool parameters, specify them in the `ext.args` entry of a process in your [custom configuration](#custom-configuration-files).

For example, if a pipeline doesn't support the `-n` parameter for `BOWTIE2_ALIGN`, add it to your custom configuration:

```groovy
process {
    withName: BOWTIE2_ALIGN {
        ext.args = "-n 0.1"
    }
}
```

When a pipeline uses a tool multiple times, specify the complete execution path of the module:

```groovy
process {
    withName: 'NFCORE_TAXPROFILER:TAXPROFILER:HOST_REMOVAL:BOWTIE2_BUILD' {
        ext.args = "-n 0.1"
    }
}
```

:::tip
Copy existing parameters from the pipeline's `conf/modules.config` file to ensure the pipeline functions correctly.
:::

:::warning
Pipeline developers provide no warranty when you update tool parameters or change `ext.args`.
:::

## Debugging

### Pipeline crashes immediately

**Problem**: Pipeline fails at the first process.

**Cause**: No container or environment profile specified.

**Solution**: Add `-profile docker`, `-profile singularity`, or `-profile conda` to your command:

```bash
nextflow run nf-core/pipeline -profile docker
```

Without a profile, Nextflow expects all tools to be manually installed on your system.

### Cluster job submission fails

**Problem**: Error message like "Failed to submit process to grid scheduler for execution" with no details.

**Solution**: Manually submit the failed job to see the actual error:

1. Navigate to the work directory shown in the error
2. Submit the command file directly to your scheduler:

   ```bash
   sbatch .command.run  # For Slurm
   qsub .command.run    # For SGE
   ```

3. The scheduler will show the real error (for example, "job memory limit exceeded")
4. Adjust your configuration based on the error message

### Invalid executor configuration

**Problem**: Error like "sbatch: error: Invalid account or account/partition combination specified".

**Common causes**:

- Missing cluster profile in `-profile` parameter
- Incorrectly specified executor in your configuration
- Resource requests exceed your cluster limits

**Solution**:

- Check if your institution has a profile in [nf-core/configs](https://github.com/nf-core/configs)
- Verify your executor name matches your scheduler (for example, `slurm` not `SLURM`)
- Confirm your `resourceLimits` match your cluster's capabilities

### Singularity bind path errors

**Problem**: "ERROR: Failed to resolve path to /home/path: No such file or directory".

**Cause**: Singularity cannot access your file system paths.

**Solution**: Add bind paths to your Nextflow configuration:

```groovy
singularity {
  enabled = true
  autoMounts = true
  runOptions = '-B /scratch -B /gpfs -B /work'
}
```

Or update your Singularity system configuration at `/etc/singularity/singularity.conf`.

### Container not updating

**Problem**: Pipeline uses old tool versions despite specifying `dev` branch.

**Cause**: Docker does not automatically update local images.

**Solution**: Manually pull the latest container:

```bash
docker pull nfcore/pipeline:dev
```

Or add to your configuration to always pull:

```groovy
docker {
  enabled = true
  runOptions = '--pull=always'
}
```

## Additional resources

For more information about configuration syntax and parameters, see:

- Links

<!-- TODO: Add links, if any -->
