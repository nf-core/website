---
title: Step-by-step guide to writing an institutional profile
subtitle: Walkthrough on what you need to set up an nf-core institutional profile
---

## Introduction

nf-core offers a centralised place to store Nextflow configuration profiles that work at an _institutional_ level.

What this means is that you can specify common nextflow pipeline configurations and options that can be shared across all users of your particular institutional cluster.

nf-core offers two level of profile sharing: institutional and pipeline-specific profiles via [nf-core/configs](https://github.com/nf-core/configs)

- Institutional configuration represent configuration options that apply to users of _all_ nf-core pipelines. These typically define settings regarding the cluster itself, such as the type of scheduler being used, maximum resource limits and so on.
- Pipeline-specific profiles that represent configuration options that apply to users of a _specific_ nf-core pipelines. These typically define common parameters all users of the pipeline would use, or customise resource requirements for particular processes of that specific pipeline

This walkthrough will guide you through making an _institutional_-level profile. It will:

- Describe commonly useful information that is worth gathering _before_ writing such a profile
- Go through each component of the profile to describe how this should be written
- Show how to test such a profile before submitting it to nf-core

## Preparation

We will first go through a list of questions to ask and what information to gather prior writing the profile. Not all questions will be relevant to your cluster, however it should cover a range of common set of cluster specifications.

### At what level should the institutional profile be designed

Various institutions have different structures in the way clusters are accessed by different groups.

This is an important consideration to make before starting writing a profile because you need to select a useful and descriptive name.

Sometimes one 'institution' will have multiple 'clusters' with the different names, which may share a common storage space. In this case it might make sense to have one nf-core configuration file for the whole institution, but with multiple _internal_ profiles representing each cluster.

Alternatively, within one 'institution', you may have groups or departments with entirely different and separated clusters with dedicated IT teams. In this case it may make sense to have separate configuration files without a common institutional 'umbrella' name.

Decide what the best structure would apply to your context.

### Do you have permission to make the profile public

nf-core institutional profiles are stored publicly on the [nf-core/configs](https://github.com/nf-core/configs/) repository.

In some cases, system administrators of the cluster at your institution may wish to keep certain aspects of the cluster private for security reasons. We therefore recommend you check with your sysadmins that you have permission to make such a profile and submit it to nf-core. We recommend you send the link to the repository with one of the configs files to show as an example of the sort of information that would be posted.

### Does your cluster use a scheduler or executor

Nextflow has integrated a range of scheduling systems to simplify running Nextflow pipelines on a range of HPCs. What this means is that it will write and sumbit submission scripts to your given scheduler or 'executor' (in Nextflow language) for you.

You can see the range of support schedulers on the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html). Check if your scheduler is there, and how it would be written (typically the name in all lower case).

You should note down if there are any special options that your cluster's scheduler requires. This maybe special ways in which memory requirements are specified, e.g.: memory per core, types of parallelisation library is required etc., the maximum number of parallel jobs that are allowed.

<details markdown="1"><summary> SGE tips</summary>

Some SGE configurations require the use of the `-l h_vmem=` flag _in addition_ to `-l virtual_free=` - the latter of which is the default used by default in Nextflow.

Furthermore, you normally need to specify the protocol to use for parallelisation (e.g. `OpenMP`, `smp`, `make`). You should find what is the recommended protocol to use.
</details>
</br>

### What container engines does your cluster offer

nf-core _highly_ recommends the use of container engines or software environment for running truly reproducible pipelines (rather than locally installed tools). This means the actual tools with compatible versions used within the pipeline are contained in a singular 'image' file.

You should find out what container engines/environments your cluster offers. For nf-core pipelines to work, you need one of any listed on the [installation](https://nf-co.re/usage/installation) page.

If you need to somehow 'load' any of the software prior use (e.g. `module load <software>` on some clusters ), you should also note that down.

Also check whether there are any special flags or parameters required for the container software itself. For example, `Singularity` sometimes may need a 'bind path' specified with `-B`.

### What are the resource maximums

All computing HPCs will have some resource limits, depending on the sizes and configuration of the nodes it offers. You should find out what is _largest_ node available to all users of a HPC and record the following information:

- maximum amount of memory (RAM) available
- maximum number of CPUs (cores) available

<details markdown="1"><summary> SLURM tips</summary>

For example, if you're using the SLURM scheduler, you can run the following command to find specifications of all nodes on your cluster

```bash
scontrol show node
```

and look for the node with largest `CPUTot` and `RealMemory` fields and node the values there.

</details>
</br>
Furthermore, you should check if there is any overall walltime (time limit) for jobs i.e. is there a maximum time a job can run for before it gets killed by the scheduler.
</br>
</br>

<details markdown="1"><summary> SLURM tips</summary>

For example, if you're using the SLURM scheduler, you can run the following command to the walltime of any queue/partitions on your cluster with

```bash
sinfo
```

and look for the TIMELIMIT column for the longest running queue that _all_ users have access to.
</details>
</br>

### Do you have to specify a queue or partition

Some clusters require your to specify which queue or partition to submit your job to (based on the resource requirements of the job). You should check if these exist, what they are called, and what limits they may have (e.g. walltimes).

You should already consider in which 'order' you would like to use these queues, to make use of the Nextflow and nf-core automated 'retry' system (i.e. if a process runs too long or uses too much memory, it can be automatically resubmitted with more time/memory)

<details markdown="1"><summary> SLURM tips</summary>

For example, if you're using the SLURM scheduler, you can run the following command to the walltime of any queues/partitions on your cluster with

```bash
sinfo
```

You should record the names of each `PARTITION` and the corresponding `TIMELIMIT` column. Note there may be additional resource limitations for each queue, so check your institutional cluster's documentation for any other limits.

</details>

## Are there any common resource or data locations on the cluster

For institutions that work on common topics (e.g. genomics), you might have a centralised local copy of the Illumina `iGenome` set of reference genomes. If it exists, you can note down the paths to this directory.

Furthermore, some institutions will specify a common place to store all container images (e.g. for Docker, Singularity etc) in a cache directory. This helps preventing multiple users having their own copies of the same pipeline container image. You should find out if, and where this location exists and note down the path.

## Should jobs be run using a local `$TMPDIR`

Some clusters require the use of a temporary _per-node_ `scratch` directory, the location of which is specified in an environmental variable `$TMPDIR`. Check if this is the case and this is set up in `$TMPDIR`

## Preparation checklist

To summarise, before you start writing the profile have you:

- Checked you know what institutional level you will make the profile at
- Checked permission to publicly post the profile
- Checked which scheduler the cluster uses
- Checked which container engine you will specify
- Checked which maximum computational resource limits exist
- Checked for any queue/partition names and specifications
- Checked for locations of any common resources or cache directories
- Checked if jobs on the cluster require use of a per-node `scratch`

## Creating the profile

Once you have all the information ready, we can begin writing the institutional profile! The tutorial will now go through each relevant Nextflow 'scope' block and explain what and when to put information in.

Before you start consider what your institutional profile should be called. The names should be as unique and recognisable as possible - however abbreviations or acronyms are fine.

For example, the name of the cluster can be used (e.g. [`phoenix`](https://github.com/nf-core/configs/blob/master/conf/phoenix.config) or [binac](https://github.com/nf-core/configs/blob/master/conf/binac.config)). In some cases an acronym of the cluster or the institute can be selected (e.g. [`SHH`](https://github.com/nf-core/configs/blob/master/conf/shh.config) or [`EVA`](https://github.com/nf-core/configs/blob/master/conf/eva.config)).

Next, make sure to make your own _fork_ of [nf-core/configs](https://github.com/nf-core/configs), and make a new branch to work on. It's normally best to call the branch the name of your proposed cluster name. Clone this to your machine, and open the folder.

### Required files

In your branch, we will need to initialise a couple of new files, and update a couple of others.

* **create** an empty file to the `conf/` directory named `<your_cluster_name>.config`
* **create** an empty file to the `docs/` directory named `<your_cluster_name>.md`
* **edit** and add your profile name to the `nfcore_custom.config` file in the top-level directory of the clone
* **edit** and add your profile name to the list of clusters on the `README.md` file in the top-level directory of the clone  under the 'Documentation'
* **edit** and add your profile name to the GitHub actions `.yaml` file (under `.github/workflows/main.yml`)

### Writing the main config file
#### params scope

In Nextflow, the `params` block of configuration files is typically used for setting pipeline-level parameters. In this case, 

define descroption, contact, nf-core max_* params, common resource locations 

#### process scope

define exectuor, queues, additional stuff e.g. penv
#### executor scope

executor limitations, e.g. maximum no. running jobs
#### <your_container> scope

container scope e.g. define which one and cache dir (+ other required options)

#### profiles{} scope

Additional profiles
### Writing the documentation file

## Testing the profile


## Additional tips and tricks
### Does your cluster often have too-full harddisk space issues

Some clusters often have

cleanup
savereference
