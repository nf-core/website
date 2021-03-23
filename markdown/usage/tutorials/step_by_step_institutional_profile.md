---
title: Step-by-step guide to writing an institutional profile
subtitle: Walkthrough on what you need to set up an nf-core institutional profile
---

## Introduction

nf-core offers a centralised place to store Nextflow configuration profiles that work at an _institutional_ level.

What this means is that you can specify common nextflow pipeline configurations and options that can be shared across all users of that particular institutional cluster.

nf-core offers two level of profile sharing: institutional and pipeline-specific profiles via [nf-core/configs](https://github.com/nf-core/configs)

- Institutional configuration represent configuration options that apply to all users of _all_ nf-core pipelines. These typically define settings regarding the cluster itself, such as the type of scheduler being used, maximum resource limits and so on.
- Pipeline-specific profiles that represent configuration options that apply to all users of a _specific_ nf-core pipelines. These typically define common parameters all users of the pipeline would use, or customise resource requirements for particular processes of that specific pipeline

This walkthrough will guide you through making an _institutional_-level profile. It will:

- Describe commonly useful information that is worth gathering _before_ writing such a profile
- Go through each component of the profile to describe how this should be written
- Show how to test such a profile before submitting it to nf-core

## Preparation

The first thing you can do is go through the following checklist

### Do you have permission to make the profile public

nf-core institutional profiles are stored publicly on the [nf-core/configs](https://github.com/nf-core/configs/) repository.

In some cases, system administrators of the cluster at your institution may which to keep certain aspects of the cluster for security reasons. We therefore recommend you check with your sysadmins that you have permission to make such a profile and submit it to nf-core. We recommend you send the link to the repository with one of the configs files to show as an example of the sort of information that would be posted.

### What container engines does your cluster offer

nf-core highly recommends the use of container engines or software environment for running truly reproducible pipelines. This means the actual tools (with correct versions) used within the pipeline are contained in a singular 'image' file.

Therefore, you need to find out what container engines/environments your cluster offers. For nf-core pipelines to work, you need one of any listed on the [installation](https://nf-co.re/usage/installation) page.

If you need to somehow 'load' any of the software prior use (e.g. `module load <software>` on some clusters ), you should also note that down.

Also check whether there are any special flags or parameters required for the container software itself. For example, `Singularity` sometimes need a 'bind path' specified with `-B`.

### Does your cluster use a scheduler or executor

Nextflow has integrated a range of scheduling systems to simplify running Nextflow pipelines on a range of HPCs. What this means is that it will write for you, and submit, submission scripts to your given scheduler or 'executor' (in Nextflow language).

You can see the range of support schedulers on the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html). You should note

You should note down if there are any special options that your cluster's schedulers requires. 

For example, for the `SGE` scheduler, some configurations require the use of the `-l h_vmem=` flag _in addition_ to `-l virtual_free=` - the latter of which is the default used by default in Nextflow.

As another example, in some cases, there may be a limit on the number of jobs that can be submitted by one user and any one time. This should also be noted down.

### What are the resource maximums

All computing HPCs will have some resource limits, depending on the sizes of the nodes it offers. You should find out what is largest node available to all users of a HPC and record the following information:

- maximum amount of memory (RAM) available
- maximum number of CPUs (cores) available

> For example, if you're using the SLURM scheduler, you can run the following command to find specifications of all nodes on your cluster
>  
>    ```bash
>    scontrol show node
>    ```
>
> and look for the node with largest `CPUTot` and `RealMemory` fields and node the values there.

Furthermore, you should check if there is any walltimes (time limit) for any jobs i.e. is there a maximum time a job can run for before it gets killed by the scheduler.

> For example, if you're using the SLURM scheduler, you can run the following command to the walltime of any queue/partitions on your cluster with
>
>    ```bash
>    sinfo
>    ```
>
> and look for the TIMELIMIT column for the longest running queue that _all_ users have access to.

### Do you have to specify a queue or partition

Some clusters require to to specify which queue or partition to submit your job to (based on the resource requirements of the job). You should check if these exist, what they are called, and what limits they may have (e.g. walltimes).

<details><summary> SLURM tips</summary>

<p>

> For example, if you're using the SLURM scheduler, you can run the following command to the walltime of any queues/partitions on your cluster with
>
>    ```bash
>    sinfo
>    ```
>
> You should record the names of each `PARTITION` and the corresponding `TIMELIMIT` column. Note there may be additional resource limitations for each queue, so check your institutional cluster's documentation for any other limits.
</p>

</details>

## Are there any common resource or data locations on the cluster

For institutions that work on common topics, they maybe have centralised data  such as the Illumina `iGenome` set of reference

## Preparation checklist

- I have permission to publicly post the profile
- I know which container engine I will specify
- I know which scheduler the cluster uses
- I know which maximum resource limits exist
- I know any queue/partition names and specifications
