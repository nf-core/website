---
title: Pipeline configuration
subtitle: How configure Nextflow to work on your system
weight: 3
---

# Introduction

One of the strongest features of Nextflow is that it can run on virtually any computational infrastructure.
It has built-in support for HPC execution schedulers such as [Slurm](https://slurm.schedmd.com/quickstart.html), [SGE](https://docs.oracle.com/cd/E19680-01/html/821-1541/ciagcgha.html#scrolltoc), [PBS](https://www.openpbs.org/), [LSF](https://www.ibm.com/support/knowledgecenter/en/SSWRJV_10.1.0/lsf_welcome/lsf_welcome.html) and more as well as cloud compute infrastructure such as [AWS Batch](https://aws.amazon.com/batch/) and [Google Cloud](https://cloud.google.com/).

Nextflow also supports container engines such as [Docker](https://www.docker.com/), [Singularity](https://sylabs.io/), [Podman](https://podman.io/), [Charliecloud](https://hpc.github.io/charliecloud/), [Shifter](https://www.nersc.gov/research-and-development/user-defined-images/), as well as the [Conda](https://docs.conda.io/en/latest/) package management system to deploy the pipelines.

In order to get nf-core pipelines to run properly on your system, you will need to install one of the aforementioned software (See [Installation](https://nf-co.re/docs/usage/installation/) page) and configure Nextflow so that it knows how best to run your analysis jobs.

# Different config locations

**You do not need to edit the pipeline code to configure nf-core pipelines.**
If you edit the pipeline defaults then you cannot update to more recent versions of the pipeline without overwriting your changes.
You also run the risk of moving away from the canonical pipeline and losing reproducibility.

There are three main types of pipeline configuration that you can use:

1. [Basic pipeline configuration profiles](#basic-configuration-profiles)
2. [Shared _nf-core/configs_ configuration profiles](#shared-nf-coreconfigs)
3. [Custom configuration files](#custom-configuration-files)

## Basic configuration profiles

Each nf-core pipeline comes with a set of "sensible defaults" for the resource requirements of each of the steps in the workflow, found in `config/base.config`.
These are always loaded and overwritten as needed by subsequent layers of configuration.

In addition to this base config, pipelines have configuration "profiles" that can be enabled with the command line flag `-profile`. Multiple profiles can be specified in a comma-separated list (e.g. `-profile test,docker`). The order of arguments is important! They are loaded in sequence, so later profiles can overwrite earlier profiles. Alternatively, you can create your own configuration profiles and supply these to the pipeline when running.

nf-core offers a range of basic profiles for configuration of container engines:

- `docker`
  - A generic configuration profile to be used with [Docker](http://docker.com/)
  - Pulls software from quay.io
- `apptainer`
  - A generic configuration profile to be used with [Singularity](https://apptainer.org/)
  - Pulls software from quay.io
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://www.nersc.gov/research-and-development/user-defined-images/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  - Pulls most software from [Bioconda](https://bioconda.github.io/)

:::note
Please only use Conda as a last resort, i.e., when it's not possible to run the pipeline with Docker or Singularity.
:::

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. **This is not recommended**.

Finally, each pipeline comes with a config profile called `test` and `test_full`. These are used for automated pipeline CI tests and will run the workflow with a minimal / full-size public dataset, respectively. They can also be used for performing test(s) run of nf-core pipeline on your infrastructure, before using your own data.

## Shared nf-core/configs

If you use a shared system with other people (such as a HPC or institutional server), it is best to use a configuration profile from [nf-core/configs](https://github.com/nf-core/configs).
These are shared config profiles loaded by all nf-core pipelines at run time.

You may find that your system already has a shared profile available here (see [https://github.com/nf-core/configs](https://github.com/nf-core/configs)).
If not, please follow the instructions in the repository README and/or the [tutorial](https://nf-co.re/docs/usage/tutorials/step_by_step_institutional_profile) to add your cluster.

## Custom configuration files

If you are the only person to be running this pipeline, you can create a local config file and use this.
Nextflow looks for these files in three locations:

1. User's home directory: `~/.nextflow/config`
2. Analysis working directory: `nextflow.config`
3. Custom path specified on the command line: `-c path/to/config` (multiple can be given)

Configuration parameters are loaded one after another and overwrite previous values.
Hardcoded pipeline defaults are first, then the user's home directory, then the work directory,
then every `-c` file in the order supplied, and finally command line `--<parameter>` options.

:::warning
For Nextflow DSL2 nf-core pipelines - parameters defined in the parameter block in `custom.config` files **WILL NOT** override defaults in `nextflow.config`! Please use `-params-file` in `yaml` or `json` format in these cases:

```json title="nf-params.json"
{
  "<parameter1_name>": 1,
  "<parameter2_name>": "<string>"
}
```

You can also generate such a JSON via each pipelines 'launch' button on the [nf-co.re website](https://nf-co.re/launch).
:::

:::warning
When tuning your pipeline configuration resources, if you want to use the `check_max()` function in your custom config, you must copy the [function](https://github.com/nf-core/tools/blob/99961bedab1518f592668727a4d692c4ddf3c336/nf_core/pipeline-template/nextflow.config#L206-L237) in the link above to the bottom of your custom config
:::

See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about configuration syntax and available parameters.

# Running Nextflow on your system

## Executor

Nextflow pipelines default to running in "local" mode unless told otherwise - that is, running jobs on the same system where Nextflow is running.
Most users will need to specify an "executor", telling Nextflow how to submit jobs to a job scheduler (_e.g._ SGE, LSF, SLURM, PBS, AWS Batch _etc._).

This can be done either within the shared configuration profiles (see [above](#shared-nf-coreconfigs)) or in custom configuration files.
For more information on how to specify an executor, please refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html#executor-page).

## Max resources

In addition to the executor, you may find that pipeline runs occasionally fail due to a particular step of the pipeline requesting more resources than you have on your system.

To avoid these failures, all nf-core pipelines [check](https://github.com/nf-core/tools/blob/99961bedab1518f592668727a4d692c4ddf3c336/nf_core/pipeline-template/nextflow.config#L206-L237) pipeline-step resource requests against parameters called `--max_cpus`, `--max_memory` and `--max_time`. These should represent the maximum possible resources of a machine or node.

These parameters only act as a _cap_, to prevent Nextflow submitting a single job requesting resources more than what is possible on your system.

:::warning
Increasing these values from the defaults will not _increase_ the resources available to the pipeline tasks! See [Tuning workflow resources](#tuning-workflow-resources) for this.
:::

Most pipelines will attempt to automatically restart jobs that fail due to lack of resources with double-requests, these caps keep those requests from getting out of hand and crashing the entire pipeline run. If a particular job exceeds the process-specific default resources and is retried, only resource requests (cpu, memory, or time) that have not yet reached the value set with `--max_<resource>` will be increased during the retry.

:::warning
Setting the `--max-<resource>` parameters do not represent the total sum of resource usage of the pipeline at a given time - only a single pipeline job!
:::

## Tuning workflow resources

The base config of nf-core pipelines defines the default resources allocated to each different step in the workflow (e.g. in a [`base.config`](https://github.com/nf-core/rnaseq/blob/master/conf/base.config) file).

These values are deliberately generous due to the wide variety of workloads done by different users.
As a result, you may find that the jobs are given more resources than they need and your system is not used as efficiently as possible.
At the other end of the scale, you may want to _increase_ the resources given to a specific task to make it run as fast as possible. You may wish to do this if you get a pipeline reporting a step failing with an `Command exit status` of e.g., `137`.

Where possible we try to get tools to make use of the resources available, for example with a tool parameter like `-p ${task.cpus}`, where `${task.cpus}` is dynamically set according to what has been made specified in the pipeline configuration files. However, this is not possible with all tools, in which case we try to make a best guess.

To tune workflow resources to better match your requirements, we can tweak these through [custom configuration files](#custom-configuration-files) or [shared nf-core/configs](#shared-nf-coreconfigs).

By default, most process resources are specified using process _labels_, for example with the following base config:

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

- The [`check_max()`](https://github.com/nf-core/tools/blob/99961bedab1518f592668727a4d692c4ddf3c336/nf_core/pipeline-template/nextflow.config#L206-L237) function applies the thresholds set in `--max_cpus`, `--max_memory` and `--max_time`.
- The `* task.attempt` means that these values are doubled if a process is automatically retried after failing with an exit code that corresponds to a lack of resources.

:::warning
If you want to use the `check_max()` function in your custom config, you must copy the function in the link above to the bottom of your custom config
:::

:::warning
You don't need to copy all of the labels into your own custom config file, only overwrite the things you wish to change
:::

If you want to give more memory to _all_ large tasks across most nf-core pipelines, would would specify in a custom config file:

```groovy
process {
  withLabel:process_high {
    memory = 200.GB
  }
}
```

You can be more specific than this by targeting a given process name instead of it's label using `withName`. You can see the process names in your console log when the pipeline is running For example:

```groovy
process {
  withName: STAR_ALIGN {
    cpus = 32
  }
}
```

In some cases, a pipeline may use a tool multiple times in the workflow. In this case you will want to specify the whole execution 'path' of the module.

```groovy
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

:::info
If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.
:::

See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#process-selectors) for more information.

Once you've written the [config file](#custom-configuration-files), and you can give this to your pipeline command with `-c`.

:::warning
Be careful with your syntax - if you set the memory to be `200` then it will get 200 _bytes_ of memory. <br/>
Memory should be in quotes with a space or without quotes using a dot: `"200 GB"` or `200.GB`. <br/>
See the Nextflow docs for [memory](https://www.nextflow.io/docs/latest/process.html#memory),
[cpus](https://www.nextflow.io/docs/latest/process.html#cpus) (int) and
[time](https://www.nextflow.io/docs/latest/process.html#time).
:::

If you think that the defaults in the pipeline are way off, please the pipeline developers know either on Slack in the channel for the pipeline or via a GitHub issue on the pipeline repository! Then we can adjust the defaults to the benefit of all pipeline users.

## Docker registries

By default, the pipelines use `quay.io` as the default docker registry for Docker and Podman images. When specifying a docker container, it will pull the image from quay.io unless you specify a full URI. For example, if the process container is:

- `biocontainers/fastqc:0.11.7--4`

By default, the image will be pulled from quay.io, resulting in a full URI of:

- `quay.io/biocontainers/fastqc:0.11.7--4`

If the `docker.registry` is specified, this will be used first. E.g. if the config value `docker.registry = 'public.ecr.aws'` is added the image will be pulled from:

- `public.ecr.aws/biocontainers/fastqc:0.11.7--4`

Alternatively, if you specify a full URI in the container specification, you can ignore the `docker.registry` setting altogether. For example, this image will be pulled exactly as specified:

- `docker.io/biocontainers/fastqc:v0.11.9_cv8`

## Updating tool versions

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of nf-core pipelines uses one container or conda environment per process which makes it much easier to maintain and update software dependencies.

If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` or `conda` definition for that process using the `withName` declaration.

For example, the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline uses a tool called Pangolin that updates an internal database of COVID-19 lineages quite frequently. It doesn't make sense to re-release the nf-core/viralrecon every time a new version of Pangolin has been released.

In this case, a user can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for the tool under `modules/nf-core/` directory of the pipeline. E.g. for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags) for Docker or [Galaxy Project](https://depot.galaxyproject.org/singularity/) for Singularity
   - Note the container version tag is identical for both container systems, but must include the 'build' ID (e.g.`--pyhdfd78af_1`)
3. Create the custom config accordingly:

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

   - For Conda (note you must check against e.g. [bioconda](https://bioconda.github.io), and this does not contain the build tag):

     ```groovy
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.1.17'
         }
     }
     ```

:::warning
It is important to note that updating containers comes with no warranty by the pipeline developers! If the update tool in the container has a major changes, this may break the pipeline.
:::

:::warning
Sometimes tool developers change how tool versions are reported between updates. Updating containers may break version reporting within the pipeline and result in missing values in MultiQC version tables.
:::

## Understanding and Modifying Tool Arguments

In some cases you may wish to understand which tool arguments or options a pipeline uses, or even update or change these for your own analyses.

You can sometimes find out what parameters are used in a tool in by checking the longer 'help' description of different pipeline parameters, e.g. by pressing the 'help' button next to [this parameter](https://nf-co.re/funcscan/1.0.1/parameters#annotation_bakta_mincontig) in [nf-core/funcscan](https://nf-co.re/funcscan).

### Finding already used arguments

However if this is not listed, there are two main places that a tool can have a tool argument specified.

Most arguments (both mandatory or optional) are defined in the `conf/modules.conf` file in the pipeline code under the `ext.args` entry. For example, you can see the default arguments used by the `SORTMERNA` step of the nf-core/rnaseq pipeline [here](https://github.com/nf-core/rnaseq/blob/6e1e448f535ccf34d11cc691bb241cfd6e60a647/conf/modules.config#LL299).

:::info
Arguments specified in `ext.args` are then inserted into the module itself via the `$args` variable in the module's bash code
:::

In some cases _some_ modules have mandatory information for a tool for it to be executed, and these normally equate to 'mandatory' arguments. You can see the argument is used in the pipeline itself looking in the `script` section given module code itself, as in the pipeline's GitHub repository under `modules/<nf-core/local>/<tool>/main.nf`.

### Customising tool arguments

If you want to modify which parameters are used by a given tool, you can do this by specifying them in the the `ext.args` entry of a process in a [custom config](#custom-configuration-files).

For example, lets say a pipeline does not yet support the `-n` parameter to a `BOWTIE_BUILD` step. You can add this in a custom config file like so:

```groovy
process {
    withName: BOWTIE2_ALIGN {
        ext.args = "-n 0.1"
    }
}
```

In some cases, a pipeline may use a tool multiple times in the workflow. In this case you will want to specify the whole execution 'path' of the module.

```groovy
process {
    withName: 'NFCORE_TAXPROFILER:TAXPROFILER:HOST_REMOVAL:BOWTIE2_BUILD' {
        ext.args = "-n 0.1"
    }
}
```

:::warning
It is recommended to copy and paste existing parameters in a pipelines `conf/modules.config` file, to ensure the pipeline can function as expected
:::

:::warning
It is important to note that updating tool parameters or changing `ext.args` comes with no warranty by the pipeline developers!
:::

# Debugging

If you have any problems configuring your profile, please see relevant sections in the [Troubleshooting documentation](https://nf-co.re/docs/usage/troubleshooting.md)
