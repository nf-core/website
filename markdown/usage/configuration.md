---
title: Nextflow configuration
subtitle: How configure Nextflow to work on your system
---

One of the strongest features of Nextflow is that it can run on virtually any computational infrastructure.
It has built-in support for HPC execution schedulers such as Slurm, SGE, PBS, LSF and more as well as cloud compute infrastructure such as AWS Batch and Google Cloud.

In order to get nf-core pipelines to run properly on your system, you will need to configure Nextflow so that it knows how best to run your analysis jobs.
Each nf-core pipeline comes with "sensible defaults" for every configuration option. These are then overwritten as needed.

Pipelines have configuration "profiles" that can be enabled with the command line flag `-profile`. Multiple profiles can be specified in a comma-separated list.
Alternatively, you can create your own configuration files and supply these to the pipeline when running.

There are three main types of pipeline configuration that you can use:

1. [Basic pipeline configuration profiles](#basic-configuration-profiles)
2. [Shared _nf-core/configs_ configuration profiles](#shared-nf-core-configs)
3. [Custom configuration files](#custom-configuration-files)

Note that most users will need to specify an "executor" - telling Nextflow how to submit jobs to a job scheduler (_e.g._ SGE, LSF, SLURM, PBS, AWS Batch _etc._).
This can be done within share configuration profiles or in custom configuration files.
For more information on how to specify an executor, please refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html#executor-page).

## Basic configuration profiles

Every nf-core pipeline comes with the config profiles `docker`, `singularity` and `conda`.
Running with these profiles will instruct the pipeline to use the relevant software packaging method.

## Shared nf-core/configs

If you use a shared system with other people, it is best to use a configuration profile from [nf-core/configs](https://github.com/nf-core/configs).
These are shared config profiles loaded by all nf-core pipelines at run time.

You may find that your system already has a shared profile available here (see [https://github.com/nf-core/configs](https://github.com/nf-core/configs)).
If not, please follow the instructions in the repository readme to add your cluster.

## Custom configuration files

If you are the only person to be running this pipeline, you can create a local config file and use this.
Nextflow looks for these files in three locations (in this order):

1. User's home directory: `~/.nextflow/config`
2. Analysis working directory: `nextflow.config`
3. Custom path specified on the command line: `-c path/to/config` (multiple can be given)

See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about configuration syntax and available parameters.
