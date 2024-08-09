---
title: System requirements
subtitle: Configure your system requirements
shortTitle: System requirements
weight: 3
---

## System requirements

nf-core pipelines are reproducible, portable, and scalable. nf-core pipelines are designed to work "out of the box". However, it is likely you will need to configure each pipeline to suit your data and system requirements.

The following sections are common considerations for nf-core pipeline users.

### Executors

Nextflow pipelines will run locally by default.
Most users will need to specify an "executor" to submit jobs to a scheduler (e.g. SGE, LSF, SLURM, PBS, and AWS Batch).

Schedulers can be specified using [shared configuration profiles](#shared-configuration-files) or [custom configuration files](#custom-parameter-and-configuration-files).
See the [Nextflow executor documentation](https://www.nextflow.io/docs/latest/executor.html#executor-page) for more information about specific schedulers.

### Max resources

Pipeline runs occasionally fail when a process requests more resources than you have available on your system.

To avoid these failures, all nf-core pipelines [check](https://github.com/nf-core/tools/blob/99961bedab1518f592668727a4d692c4ddf3c336/nf_core/pipeline-template/nextflow.config#L206-L237) pipeline-step resource requests against parameters called `--max_cpus`, `--max_memory` and `--max_time`. These parameters can be set as the maximum possible resources of a machine or node and act as a _cap_ that prevents Nextflow submitting a single job requesting resources more than what is possible.

:::warning
Increasing these values from the defaults will not _increase_ the resources available to the pipeline tasks! See [Tuning resources](#tuning-resources) for this.
:::

Most pipelines will attempt to automatically restart jobs that fail due to lack of resources with double-requests. These `--max_<resource>` caps keep requests from getting out of hand and crashing the entire pipeline run. If a particular job exceeds the process-specific default resources and is retried, only resource requests (cpu, memory, or time) that have not yet reached the value set with `--max_<resource>` will be increased during the retry.

:::warning
The `--max_<resource>` parameters do not represent the total sum of resource usage of the pipeline at a given time - only a single pipeline job!
:::

### Tuning resources

The base config of nf-core pipelines defines the default resources allocated to each different step in the workflow (e.g., in a [`base.config`](https://github.com/nf-core/rnaseq/blob/master/conf/base.config) file).

These values are deliberately generous due to the wide variety of workloads done by different users.
As a result, you may find that the jobs are given more resources than they need and your system is not used efficiently.
At the other end of the scale, you may want to increase the resources given to a specific task to make it run faster.
You may wish to do this if you get a pipeline reporting a step failing with an `Command exit status`, such as `137`.

Where possible tools are tuned to make use of the resources available, for example with a tool parameter (e.g., `-p ${task.cpus}`, where `${task.cpus}` is dynamically set according to what has been specified in the pipeline configuration files). However, this is not possible with all tools, in which case required resources are estimated for an average user.

Workflow resources can modified through [shared configuration profiles](#shared-configuration-files) or [custom configuration files](#custom-parameter-and-configuration-files) to better match your requirements.

Most process resources are specified using process _labels_ by default. For example:

```groovy
process {
  withLabel:process_low {
    cpus = { check_max( 2 * task.attempt, 'cpus' ) }
    memory = { check_max( 14.GB * task.attempt, 'memory' ) }
    time = { check_max( 6.h * task.attempt, 'time' ) }
  }
  withLabel:process_medium {
    cpus = { check_max( 6 * task.attempt, 'cpus' ) }
    memory = { check_max( 42.GB * task.attempt, 'memory' ) }
    time = { check_max( 8.h * task.attempt, 'time' ) }
  }
  withLabel:process_high {
    cpus = { check_max( 12 * task.attempt, 'cpus' ) }
    memory = { check_max( 84.GB * task.attempt, 'memory' ) }
    time = { check_max( 10.h * task.attempt, 'time' ) }
  }
}
```

The [`check_max()`](https://github.com/nf-core/tools/blob/99961bedab1518f592668727a4d692c4ddf3c336/nf_core/pipeline-template/nextflow.config#L206-L237) function applies the thresholds set in `--max_cpus`, `--max_memory` and `--max_time`. The `* task.attempt` means that these values are doubled and automatically retried after failing with an exit code that corresponds to a lack of resources.

:::warning
If you want to use the `check_max()` function in a custom configuration file, you must copy the [check_max function](https://github.com/nf-core/tools/blob/99961bedab1518f592668727a4d692c4ddf3c336/nf_core/pipeline-template/nextflow.config#L206-L237) to the bottom of your custom config
:::

To modify the memory to all processes with the `process_high` label you would use the `withLabel` process selector. For example:

```groovy
process {
  withLabel:process_high {
    memory = 200.GB
  }
}
```

To modify the memory of a specific process with the name `STAR_ALIGN` you would use the `withName` process selector. For example:

```groovy
process {
  withName: STAR_ALIGN {
    memory = 200.GB
  }
}
```

If a pipeline uses a tool multiple times you will need to specify the whole execution path of the module. For example:

```groovy
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#process-selectors) for more information.

:::warning
Be careful with your syntax - if you set the memory to be `200` then it will get 200 _bytes_ of memory. <br/>
Memory should be in quotes with a space or without quotes using a dot: `"200 GB"` or `200.GB`. <br/>
See the Nextflow docs for [memory](https://www.nextflow.io/docs/latest/process.html#memory),
[cpus](https://www.nextflow.io/docs/latest/process.html#cpus) (int) and
[time](https://www.nextflow.io/docs/latest/process.html#time) for more information.
:::

:::note
If you think that the default resources for a pipeline are drastically too high or low please contact the developers know either on Slack in the channel for the pipeline or via a GitHub issue on the pipeline repository.
:::
