---
title: Pipeline configuration
subtitle: How configure Nextflow to work on your system
---

- [Introduction](#introduction)
- [Different config locations](#different-config-locations)
  - [Basic configuration profiles](#basic-configuration-profiles)
  - [Shared nf-core/configs](#shared-nf-core-configs)
  - [Custom configuration files](#custom-configuration-files)
- [Running Nextflow on your system](#running-nextflow-on-your-system)
  - [Executor](#executor)
  - [Max resources](#max-resources)
  - [Tuning workflow resources](#tuning-workflow-resources)

# Introduction

One of the strongest features of Nextflow is that it can run on virtually any computational infrastructure.
It has built-in support for HPC execution schedulers such as Slurm, SGE, PBS, LSF and more as well as cloud compute infrastructure such as AWS Batch and Google Cloud.

In order to get nf-core pipelines to run properly on your system, you will need to configure Nextflow so that it knows how best to run your analysis jobs.

# Different config locations

**You do not need to edit the pipeline code to configure nf-core pipelines.**
If you edit the pipeline defaults then you cannot update to more recent versions of the pipeline without overwriting your changes.
You also run the risk of moving away from the canonical pipeline and losing reproducibility.

Pipelines have configuration "profiles" that can be enabled with the command line flag `-profile`.
Multiple profiles can be specified in a comma-separated list.
Alternatively, you can create your own configuration files and supply these to the pipeline when running.

There are three main types of pipeline configuration that you can use:

1. [Basic pipeline configuration profiles](#basic-configuration-profiles)
2. [Shared _nf-core/configs_ configuration profiles](#shared-nf-coreconfigs)
3. [Custom configuration files](#custom-configuration-files)

## Basic configuration profiles

Each nf-core pipeline comes with a set of "sensible defaults" for the resource requirements of each of the steps in the workflow, found in `config/base.config`.
These are always loaded and overwritten as needed by subsequent layers of configuration.

In addition to this base config, pipelines come with basic config profiles to enable software packaging methods: `docker`, `singularity`, `podman` and `conda`.

Finally, each pipeline comes with a config profile called `test` and `test_full`. These are used for automated pipeline CI tests and will run the workflow with a minimal / full-size public dataset, respectively.

## Shared nf-core/configs

If you use a shared system with other people, it is best to use a configuration profile from [nf-core/configs](https://github.com/nf-core/configs).
These are shared config profiles loaded by all nf-core pipelines at run time.

You may find that your system already has a shared profile available here (see [https://github.com/nf-core/configs](https://github.com/nf-core/configs)).
If not, please follow the instructions in the repository readme to add your cluster.

## Custom configuration files

If you are the only person to be running this pipeline, you can create a local config file and use this.
Nextflow looks for these files in three locations:

1. User's home directory: `~/.nextflow/config`
2. Analysis working directory: `nextflow.config`
3. Custom path specified on the command line: `-c path/to/config` (multiple can be given)

Configuration parameters are loaded one after another and overwrite previous values.
Hardcoded pipeline defaults are first, then the user's home directory, then the work directory,
then every `-c` file in the order supplied, and finally command line `--parameter` options.

See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about configuration syntax and available parameters.

# Running Nextflow on your system

## Executor

Nextflow pipelines default to running in "local" mode unless told otherwise - that is, running jobs on the same system where Nextflow is running.
Most users will need to specify an "executor", telling Nextflow how to submit jobs to a job scheduler (_e.g._ SGE, LSF, SLURM, PBS, AWS Batch _etc._).

This can be done within share configuration profiles or in custom configuration files.
For more information on how to specify an executor, please refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html#executor-page).

## Max resources

In addition to the executor, you may find that pipeline runs occasionally fail due to a step requesting more resources than you have on your system.
To avoid these failures, all nf-core pipelines check parameters called `--max_cpus`, `--max_memory` and `--max_time`.

**:warning: Increasing these values from the defaults will not increase the resources available to the pipeline tasks!**

These parameters only act as a _cap_, to prevent Nextflow from going over what is possible on your system.
Most pipelines will attempt to automatically restart jobs that fail due to lack of resources with double-requests,
these caps keep those requests from getting out of hand and crashing the entire pipeline run.

## Tuning workflow resources

The base config of nf-core pipelines defines the default resources allocated to each different step in the workflow.
These values are deliberately generous due to the wide variety of workloads done by different users.
As a result, you may find that the jobs are given more resources than they need and your system is not used as efficiently as possible.
At the other end of the scale, you may want to _increase_ the resources given to a specific task to make it run as fast as possible.

Where possible we try to get tools to make use of the resources available, for example with a parameter like `-p ${task.cpus}` (where `${task.cpus}` is dynamically set according to what has been made available in the pipeline config). However, this is not possible with all tools, in which case we try to make a best guess.

Most process resources are specified using process _labels_, for example with the following base config:

```nextflow
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

> The `check_max()` function applies the thresholds set in `--max_cpus`, `--max_memory` and `--max_time`.

> The `* task.attempt` means that these values are doubled if a process is automatically retried after failing with an exit code that corresponds to a lack of resources.

**:warning: You don't need to copy all of this into your own custom config file, only overwrite the things you wish to change**

For example, if you want to give more memory to all large tasks across most nf-core pipelines, the following config could work:

```nextflow
process {
  withLabel:process_high {
    memory = 200.GB
  }
}
```

You can be more specific than this by targetting a given process name instead of it's label using `withName`. For example:

```nextflow
process {
  withName:bwa_align {
    cpus = 32
  }
}
```

See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#process-selectors) for more information.

> :warning: Be careful with your syntax - if you set the memory to be `200` then it will get 200 _bytes_ of memory.
> Memory should be in quotes with a space or without quotes using a dot: `"200 GB"` or `200.GB`.
> See the Nextflow docs for [memory](https://www.nextflow.io/docs/latest/process.html#memory),
> [cpus](https://www.nextflow.io/docs/latest/process.html#cpus) (int) and
> [time](https://www.nextflow.io/docs/latest/process.html#time).

Usually you would save configuration such as this to a config file [as described above](#custom-configuration-files).

If you think that the defaults in the pipeline are way off, please let us know! Then we can adjust the defaults to the benefit of all pipeline users.
