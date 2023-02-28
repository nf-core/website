---
title: Step-by-step guide to writing an institutional profile
subtitle: Walkthrough on what you need to set up an nf-core institutional profile
---

## Introduction

🏫 nf-core offers a centralised place to store Nextflow configuration profiles that work at an _institutional_ level.

👨‍👩‍👧‍ What this means is that you can specify common Nextflow pipeline configurations and options that can be shared across all users of your particular institutional cluster.

nf-core offers two level of profile sharing: global institutional and pipeline institutional profiles via [nf-core/configs](https://github.com/nf-core/configs)

- **Global** institutional profiles represent configuration options that apply to users of _all_ nf-core pipelines. These typically define settings regarding the cluster itself, such as the type of scheduler being used, maximum resource limits and so on.
- **Pipeline** institutional profiles that represent configuration options that apply to users of a _specific_ nf-core pipelines. These typically define common parameters all users of the pipeline would use, or customise resource requirements for particular processes of that specific pipeline

This walkthrough will guide you through making an _institutional_-level profile. It will:

- 🔍 Describe commonly useful information that is worth gathering _before_ writing such a profile
- ✍ Go through each component of the profile to describe how this should be written
- ⌨️ Show how to test such a profile before submitting it to nf-core

## Preparation

We will first go through a list of questions to ask and what information to gather prior writing the profile. Not all questions will be relevant to your cluster, however it should cover a range of common set of cluster specifications.

### Does a institutional profile already exist

First you should check nf-core/configs that a global institutional profile already exists!

If one exists, try it out and you can jump straight into running nf-core pipelines 😎

If you see any issues with it, or it's not working for you - contact the person who originally made it.

If no profile already exists, continue with this walkthrough!

### At what level should the global institutional profile be designed

Various institutions have different structures in the way clusters are accessed by different groups.

This is an important consideration to make before starting writing a profile because you need to select a useful and descriptive name.

Sometimes one institution will have multiple clusters with the different names, which may share a common storage space. In this case it might make sense to have one 'umbrella' nf-core global profile for the whole institution, but with multiple _internal_ profiles representing each cluster.

Alternatively, within one institution, you may have groups or departments with entirely different and separated clusters with dedicated IT teams. In this case it may make sense to have separate institutional profiles (and therefore not under common institutional 'umbrella' name).

Decide what the best structure would apply to your context.

### Do you have permission to make the profile public

nf-core global institutional profiles are stored publicly on the [nf-core/configs](https://github.com/nf-core/configs/) repository.

In some cases, system administrators of the cluster at your institution may wish to keep certain aspects of the cluster private for security reasons.

We therefore recommend you check with your sysadmins that you have permission to make such a profile before submitting it to nf-core. We recommend you send the link to the repository with one of the already existing profiles to show as an example of the sort of information that would be posted.

### Does your cluster use a scheduler or executor

Nextflow has integrated a range of scheduling systems to simplify running Nextflow pipelines on a range of HPCs.

What this means is that it will write and submit submission scripts to your given scheduler or 'executor' (in Nextflow language) for you 🤖.

You can see the range of support schedulers on the [Nextflow documentation](https://www.Nextflow.io/docs/latest/executor.html). Check if your scheduler is there, and how it would be written (typically the name in all lower case).

You should note down if there are any special options that your cluster's scheduler requires. This maybe special ways in which memory requirements are specified, e.g.: memory per core, types of parallelisation library is required etc., the maximum number of parallel jobs that are allowed ⇉.

> 🌞 Some SGE configurations require the use of the `-l h_vmem=` flag _in addition_ to `-l virtual_free=` - the latter of which is used by default in Nextflow. Furthermore, you normally need to specify the protocol to use for parallelisation (e.g. `OpenMP`, `smp`, `make`). You should find what is the recommended protocol to use.

### What container engines does your cluster offer

nf-core _highly_ recommends the use of container engines or software environment for running truly reproducible pipelines (rather than locally installed tools). This means the actual tools with compatible versions used within the pipeline are contained in a singular 'image' file.

You should find out what container engines/environments your cluster offers. For nf-core pipelines to work, you need one of any listed on the [installation](https://nf-co.re/docs/usage/installation) page.

If you need to somehow 'load' any of the software prior use (e.g. `module load <software>` on some clusters), you should also note that down.

Also check whether there are any special flags or parameters required for the container software itself. For example, `Singularity` sometimes may need a 'bind path' specified with `-B`.

### What are the resource maximums

All computing HPCs will have some resource limits, depending on the sizes and configuration of the nodes it offers. You should find out what is _largest_ node available to all users of a HPC and record the following information:

- maximum amount of memory (RAM) available
- maximum number of CPUs (cores) available

> 🐛 For example, if you're using the SLURM scheduler, you can run the following command to find specifications of all nodes on your cluster
>
> ```bash
> scontrol show node
> ```
>
> and look for the node with largest `CPUTot` and `RealMemory` fields and node the values there.

Furthermore, you should check if there is any overall walltime (time limit) for jobs i.e. is there a maximum time a job can run for before it gets killed by the scheduler.

> 🐛 For example, if you're using the SLURM scheduler, you can run the following command to the walltime of any queue/partitions on your cluster with
>
> ```bash
> sinfo
> ```
>
> and look for the TIMELIMIT column for the longest running queue that _all_ users have access to.

### Do you have to specify a queue or partition

Some clusters require you to specify which queue or partition to submit your job to (based on the resource requirements of the job). You should check if these exist, what they are called, and what limits they may have (e.g. walltimes).

You should already consider in which 'order' you would like to use these queues, to make use of the Nextflow and nf-core automated 'retry' system (i.e. if a process runs too long or uses too much memory, it can be automatically resubmitted with more time/memory)

> 🐛 For example, if you're using the SLURM scheduler, you can run the following command to the walltime of any queues/partitions on your cluster with
>
> ```bash
> sinfo
> ```
>
> You should record the names of each `PARTITION` and the corresponding `TIMELIMIT` column. Note there may be additional resource limitations for each queue, so check your institutional cluster's documentation for any other limits.

### Are there any common resource or data locations on the cluster

For institutions that work on common topics (e.g. genomics), you might have a centralised local copy of the AWS [`iGenomes`](https://ewels.github.io/AWS-iGenomes/) set of reference genomes. If it exists, you can note down the paths to this directory.

Furthermore, some institutions will specify a common place to store all container images (e.g. for Conda environments, Singularity images etc) in a cache directory. This helps preventing multiple users having their own copies of the same pipeline container image. You should find out if, and where this location exists and note down the path.

### Should jobs be run using a local `$TMPDIR`

Some clusters require the use of a temporary _per-node_ `scratch` directory, the location of which is specified in an environmental variable `$TMPDIR`. Check if this is the case and this is set up in `$TMPDIR`

### Preparation checklist

To summarise, before you start writing the profile have you checked:

- 🏫 At what institutional level you will make the global institutional profile
- 🤫 The permission to publicly post the global institutional profile
- 🔀 Which scheduler the cluster uses
- 📦 Which container engine you will specify
- 📊 Which maximum computational resource limits exist
- 📛 Queue/partition names and specifications
- 📁 Locations of any common resources or cache directories
- 📁 If jobs on the cluster require use of a per-node `scratch`

## Creating the profile

Once you have all the information ready, we can begin writing the global institutional profile! The walkthrough will now go through each relevant Nextflow 'scope' block and explain what and when to put information in.

Before you start consider what your global institutional profile should be called. The names should be as unique and recognisable as possible - however abbreviations or acronyms are fine.

For example, the name of the cluster can be used (e.g. [`phoenix`](https://github.com/nf-core/configs/blob/master/conf/phoenix.config) or [`binac`](https://github.com/nf-core/configs/blob/master/conf/binac.config)). In some cases an acronym of the cluster or the institute can be selected (e.g. [`EVA`](https://github.com/nf-core/configs/blob/master/conf/eva.config)).

Next, we suggest to make your own _fork_ of [nf-core/configs](https://github.com/nf-core/configs) (although this is not required), and make a new branch to work on. It's normally best to call the branch the name of your proposed cluster name. Clone this to your machine, and open the folder.

### Required files

In your branch, we will need to initialise a couple of new files, and update a couple of others.

- **create** an empty file to the `conf/` directory named `<your_cluster_name>.config`
- **create** an empty file to the `docs/` directory named `<your_cluster_name>.md`
- **edit** and add your profile name to the `nfcore_custom.config` file in the top-level directory of the clone
- **edit** and add your profile name to the list of clusters on the `README.md` file in the top-level directory of the clone under the 'Documentation'
- **edit** and add your profile name to the GitHub actions `.yaml` file (under `.github/workflows/main.yml`)

### Writing the main global institutional profile file

First we will edit the main profile file under `conf/<your_cluster_name>.config`.

#### params scope

In Nextflow, the `params` block of configuration files is typically used for setting pipeline-level parameters. In the case of global institutional profiles we will very likely not specify pipeline parameters here, but rather add some useful nf-core specific parameters that apply to all pipelines. See the nf-core/configs README for more information how to define _pipeline_ institutional profiles.

The most useful first step for testing a new nf-core global institutional profile is to add to the params scope the `config_profile_*` series of params.

These are nf-core specific parameters that are displayed in the header summary of each run, describing what the profile is and who maintains it. Therefore, you can use this when testing the profile to check the profile was actually loaded.

In the `conf/<your_cluster_name>.config` file, add to a params scope something like the following:

```nextflow
params {
  config_profile_description = '<your_cluster_name> cluster profile provided by nf-core/configs.'
  config_profile_contact = '<your_name> (<your_github_handle>)'
  config_profile_url = 'https://<institutional_url>.com'
}
```

Note that for the `config_profile_contact`, it is best to indicate a specific person. This will typically be someone who wrote the config (via their name & github handle) or whoever will maintain it at the institution (e.g. email of IT Department, Institution X), i.e. someone who can be contacted if there are questions or problems and how to contact them.

Next, in the same scope, we can also specify the `max_*` series of params.

These are used by nf-core pipelines to limit automatic resubmission of resource-related failed jobs to ensure submitted retries do not exceed the maximum available on your cluster. These values should be the ones you found for the largest node of your cluster (i.e., the largest node a user's job can be submitted to). For example:

```nextflow
params {
  config_profile_description = '<your_cluster_name> cluster profile provided by nf-core/configs.'
  config_profile_contact = '<your_name> (<your_github_handle>)'
  config_profile_url = 'https://<institutional_url>.com'
  max_memory = 2.TB
  max_cpus = 128
  max_time = 720.h
}
```

Finally, if you have a common resource directory for the AWS `iGenomes` collection of reference genomes, this can can also go in the `params` scope.

```nextflow
params {
  config_profile_description = '<your_cluster_name> cluster profile provided by nf-core/configs.'
  config_profile_contact = '<your_name> (<your_github_handle>)'
  config_profile_url = 'https://<institutional_url>.com'
  max_memory = 2.TB
  max_cpus = 128
  max_time = 720.h
  igenomes_base = '/<path>/<to>/igenomes/'
}
```

#### process scope

Next, we can use the `process` scope to define which scheduler to use and associated options. Any option specified in this scope means that all processes in a pipeline will use the settings defined here.

Normally, you only need to specify which scheduler you use. For example, if using SLURM 🐛:

```nextflow
params {
  config_profile_description = '<your_cluster_name> cluster profile provided by nf-core/configs.'
  config_profile_contact = '<your_name> (<your_github_handle>)'
  config_profile_url = 'https://<institutional_url>.com'
  max_memory = 2.TB
  max_cpus = 128
  max_time = 720.h
  igenomes_base = '/<path>/<to>/igenomes/'
}

process {
  executor = 'slurm'
}

```

If you need to specify more cluster-specific information regarding your cluster, this can also go here.

For example, you can specify which queue to use. If you only have a single queue, this is as simple as:

```nextflow
params {
  config_profile_description = '<your_cluster_name> cluster profile provided by nf-core/configs.'
  config_profile_contact = '<your_name> (<your_github_handle>)'
  config_profile_url = 'https://<institutional_url>.com'
  max_memory = 2.TB
  max_cpus = 128
  max_time = 720.h
  igenomes_base = '/<path>/<to>/igenomes/'
}

process {
  executor = 'slurm'
  queue = 'all'
}

```

Alternatively, if you must specify different queues based on various resource-related specifications, this can be added with a Groovy expression ✌.
Lets say you have three queues that act as priority queues based on maximum runtime of jobs (short - 2 hour walltime, medium - 24 hour walltime, long - no walltime), you can specify this as follows:

```nextflow
params {
  config_profile_description = '<your_cluster_name> cluster profile provided by nf-core/configs.'
  config_profile_contact = '<your_name> (<your_github_handle>)'
  config_profile_url = 'https://<institutional_url>.com'
  max_memory = 2.TB
  max_cpus = 128
  max_time = 720.h
  igenomes_base = '/<path>/<to>/igenomes/'
}

process {
  executor = 'slurm'
  queue = { task.time <= 2.h ? 'short' : task.time <= 24.h ? 'medium': 'long' }
}

```

Alternatively, if it's based on CPU requirements of nodes:

```nextflow
params {
  config_profile_description = '<your_cluster_name> cluster profile provided by nf-core/configs.'
  config_profile_contact = '<your_name> (<your_github_handle>)'
  config_profile_url = 'https://<institutional_url>.com'
  max_memory = 2.TB
  max_cpus = 128
  max_time = 720.h
  igenomes_base = '/<path>/<to>/igenomes/'
}

process {
  executor = 'slurm'
  queue = { task.cpus > 24 ? 'big' : 'small' }
}

```

As mentioned above, Nextflow has a clever automated retry system where if a particular submission exits with certain resource-limit reached exit codes, a resubmission will be made with greater resource requests.

To ensure you can use this but also exploit the additional nf-core checks that prevent jobs from requesting more than available on a given node (and resulting the pipeline stalling because jobs are stuck in a queue), you should also specify a maximum number of retries with the Nextflow `maxRetries` directives.

```nextflow
params {
  config_profile_description = '<your_cluster_name> cluster profile provided by nf-core/configs.'
  config_profile_contact = '<your_name> (<your_github_handle>)'
  config_profile_url = 'https://<institutional_url>.com'
  max_memory = 2.TB
  max_cpus = 128
  max_time = 720.h
  igenomes_base = '/<path>/<to>/igenomes/'
}

process {
  executor = 'slurm'
  queue = { task.cpus > 24 ? 'big' : 'small' }
  maxRetries = 2
}

```

In this case, after the initial submission of a job, on resource-related failures Nextflow will retry just 2 more times before the pipeline as a whole will fail.

If you normally need to specify additional 'non-standard' options in the headers of scheduler batch scripts (e.g. `sbatch` for SLURM 🐛 or `qsub` for SGE 🌞), you can specify these with `clusterOptions`. Anything specified in the `clusterOptions` directive will be added in the header of the Nextflow-generated batch script for you (you can see these in the `.command.run` file in each job's `work/<hash>` directory in a Nextflow run).

> 🌞 For example, for some SGE clusters, memory requests are specified with the `h_vmem` variable, rather than the Nextflow default `virtual_free`.
>
> ```nextflow
> params {
>  config_profile_description = '<your_cluster_name> cluster profile provided by nf-core/configs.'
>   config_profile_contact = '<your_name> (<your_github_handle>)'
>   config_profile_url = 'https://<institutional_url>.com'
>   max_memory = 2.TB
>   max_cpus = 128
>   max_time = 720.h
>   igenomes_base = '/<path>/<to>/igenomes/'
> }
>
> process {
>   executor = 'sge'
>   queue = { task.cpus > 24 ? 'big' : 'small' }
>   maxRetries = 2
>   clusterOptions = { "-l h_vmem=${task.memory.toGiga()}G" }
> }
> ```

Where the pipeline-defined memory specification of the each job is inserted into the batch script header using the Nextflow `${task.memory}` variable.

Another commonly used directive is the `beforeScript` directive, which allows you to run a custom unix command _prior_ to running a pipeline's command of a particular job. This is often used when a UNIX software module needs to be loaded on the node the job is sent to by the scheduler.

For example, you may need to explicitly load the `singularity` container software module, which can be specified like so:

```nextflow
params {
  config_profile_description = '<your_cluster_name> cluster profile provided by nf-core/configs.'
  config_profile_contact = '<your_name> (<your_github_handle>)'
  config_profile_url = 'https://<institutional_url>.com'
  max_memory = 2.TB
  max_cpus = 128
  max_time = 720.h
  igenomes_base = '/<path>/<to>/igenomes/'
}

process {
  executor = 'sge'
  queue = { task.cpus > 24 ? 'big' : 'small' }
  maxRetries = 2
  clusterOptions = { "-l h_vmem=${task.memory.toGiga()}G" }
  beforeScript = 'module load singularity'
}
```

For a full list of `process` directives, please see the [Nextflow documentation](https://www.nextflow.io/docs/latest/process.html#process-directives)

#### executor scope

The executor scope allows the use of further executor _specific_ options that are inbuilt into Nextflow. You should check the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#scope-executor) to see what options are available for your respective scheduler.

One sometimes useful option for smaller clusters with less sophisticated fair-use management is the `queueSize` directive. This allows you to specify the maximum number of jobs of a given Nextflow run can submit in parallel at any one time. So to prevent Nextflow from swamping a (small) cluster thousands of jobs at once and blocking the cluster for other users 😱, you can limit this as follows:

```nextflow
params {
  config_profile_description = '<your_cluster_name> cluster profile provided by nf-core/configs.'
  config_profile_contact = '<your_name> (<your_github_handle>)'
  config_profile_url = 'https://<institutional_url>.com'
  max_memory = 2.TB
  max_cpus = 128
  max_time = 720.h
  igenomes_base = '/<path>/<to>/igenomes/'
}

process {
  executor = 'sge'
  queue = { task.cpus > 24 ? 'big' : 'small' }
  maxRetries = 2
  clusterOptions = { "-l h_vmem=${task.memory.toGiga()}G" }
  beforeScript = 'module load singularity'
}

executor {
  queueSize = 8
}
```

Where a given Nextflow run can only have 8 submitted jobs at once (and will wait to one job is completed before submitting the next one).

A similar directive is the `submitRateLimit` which specifies how many jobs can be specified in a given time frame (as some clusters may penalise you for over-submitting):

```nextflow
params {
  config_profile_description = '<your_cluster_name> cluster profile provided by nf-core/configs.'
  config_profile_contact = '<your_name> (<your_github_handle>)'
  config_profile_url = 'https://<institutional_url>.com'
  max_memory = 2.TB
  max_cpus = 128
  max_time = 720.h
  igenomes_base = '/<path>/<to>/igenomes/'
}

process {
  executor = 'sge'
  queue = { task.cpus > 24 ? 'big' : 'small' }
  maxRetries = 2
  clusterOptions = { "-l h_vmem=${task.memory.toGiga()}G" }
  beforeScript = 'module load singularity'
}

executor {
  queueSize = 8
  submitRateLimit = '10 sec'
}
```

Where a maximum of 10 jobs can be submitted per second.

#### container scopes

You can specify institutional cluster specific options for the container engine system (e.g., Singularity, Docker, Podman, etc.) that you will use with a variety of different container scopes . There is generally one scope per container engine (and `conda`), and more configuration information on each engine offered by Nextflow can be seen in the [Nextflow Documentation](https://www.nextflow.io/docs/latest/config.html#config-scopes).

To use the example of Singularity, you can use the singularity scope to specify that it should be used, that it should be automatically mounted and also (where valid) where a common cache directory of images resides:

```nextflow
params {
  config_profile_description = '<your_cluster_name> cluster profile provided by nf-core/configs.'
  config_profile_contact = '<your_name> (<your_github_handle>)'
  config_profile_url = 'https://<institutional_url>.com'
  max_memory = 2.TB
  max_cpus = 128
  max_time = 720.h
  igenomes_base = '/<path>/<to>/igenomes/'
}

process {
  executor = 'sge'
  queue = { task.cpus > 24 ? 'big' : 'small' }
  maxRetries = 2
  clusterOptions = { "-l h_vmem=${task.memory.toGiga()}G" }
  beforeScript = 'module load singularity'
}

executor {
  queueSize = 8
  submitRateLimit = '10 sec'
}

singularity {
  enabled = true
  autoMounts = true
  cacheDir = '/<path>/<to>/<your>/<image_cache>'
}
```

Each container engine or software environment may have different options, so be sure to check the Nextflow documentation what options you may have.

#### profiles scope

In some cases, you may want to define multiple contexts that have different specifications - For example, your institution may have has two distinct clusters or set of nodes with slightly different specifications, while still sharing the same storage.

One concise way of specifying these is via the `profiles` scope.

These essentially act as nested or internal profiles, where you can extend or modify the base parameters set in the main config.

Using our example above, maybe our institution has two clusters named red and blue, one that has larger nodes than the other. What we can do here is specify the `max_*` series of parameters in cluster-specific internal profiles.

```nextflow
process {
  executor = 'sge'
  queue = { task.cpus > 24 ? 'big' : 'small' }
  maxRetries = 2
  clusterOptions = { "-l h_vmem=${task.memory.toGiga()}G" }
  beforeScript = 'module load singularity'
}

executor {
  queueSize = 8
  submitRateLimit = '10 sec'
}

singularity {
  enabled = true
  autoMounts = true
  cacheDir = '/<path>/<to>/<your>/<image_cache>'
}

profiles {
  red {
    params {
      config_profile_description = '<your_institution_name> 'red' cluster cluster profile provided by nf-core/configs.'
      config_profile_contact = '<your_name> (<your_github_handle>)'
      config_profile_url = 'https://<institutional_url>.com'
      max_memory = 2.TB
      max_cpus = 128
      max_time = 720.h
      igenomes_base = '/<path>/<to>/igenomes/'
    }
  }

  blue {
    params {
      config_profile_description = '<your_institution_name> 'blue' cluster profile provided by nf-core/configs.'
      config_profile_contact = '<your_name> (<your_github_handle>)'
      config_profile_url = 'https://<institutional_url>.com`'
      max_memory = 256.GB
      max_cpus = 64
      max_time = 24.h
      igenomes_base = '/<path>/<to>/igenomes/'
    }
  }
}
```

You can see here we have moved the `params` block into each of the _internal profiles_, and updated the `config_profile_description` and `max_*` parameters accordingly.

> :warning: Important: you should **not** define scopes both in the global profile AND in the internal profile. Internal profiles do _not_ inherit directives/settings defined in scopes in the base config, so anything defined in the base global profile file will be _ignored_ in the internal profile. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles) for more information.

### Writing the documentation file

Once you've written your configuration file and saved it, you should briefly describe your global institutional profile in your documentation markdown file.

You can look at some of the other configs to see what this describes, but typically it will say where exactly the clusters are based, how to user can use the global institutional profile (once eventually merged into the main nf-core/configs repository), and any other particular information that would be relevant for users of the profile (such as what internal profiles are offered, or if users have to do any other personal set-up of certain variables).

## Testing the global institutional profile

Once all of your files are prepared, it's time to test your global institutional profile. All nf-core pipelines come with a special parameter called `--custom_config_base`. This allows you to override the pipeline run from loading the nf-core/configs from the nf-core repository, and you can instead get it to look in your fork.

> Note: At this stage don't get disheartened if you hit errors and problems along the way. This often requires a lot of trial and error, as you learn about the particulars of your particular cluster. All are different, and all require TLC to get everything optimal!

To test your config, pick your pipeline of choice, and run the given pipeline's integrated mini-test profile but using your config. For example:

```bash
nextflow run nf-core/eager \
-profile <your_cluster_name>,test_tsv \
--custom_config_base 'https://raw.githubusercontent.com/<your_github_user>/configs/<your_branch>'
```

If you also wish to test a internal profile from your `profiles` scope, include that in the `-profile` flag. For example:

```bash
nextflow run nf-core/eager \
-profile <your_cluster_name>,<you_profile_name>,test_tsv \
--custom_config_base 'https://raw.githubusercontent.com/<your_github_user>/configs/<your_branch>'
```

To know if either the global or internal profiles is working you can check for things like the following:

- Does the pipeline run summary displayed in your terminal at the beginning of the run display the right `config_profile_*` parameters?
  - Note in some cases the _description_ may be of the test profile. However the URL should always be of your global institutional profile: e.g.
- Does the pipeline run summary display the correct `containerEngine`?
- Does your scheduler report jobs have been submitted to it from your pipeline?
  - You can check this with `squeue` in SLURM or `qstat` in SGE, for example
  - You should see job names beginning with `nf_` in your submission log

> 💡 Tip: If you can't see the run summary header due to a higher number process status bars, you can run `export NXF_ANSI_LOG=false` before running the test command to use a more condense report without the fancy status information.

## Make a PR into nf-core configs

Once your testing works without any errors, it's time to make your global institutional profile official!

Simply make a PR into the nf-core/configs repository, and request a review.

🎉 Once approved and merged, let your colleagues know that any user can immediately use the global institutional profile for all nf-core pipelines with just the `-profile` flag!

```bash
nextflow run nf-core/<pipeline> \
-profile <your_cluster_name> \
[...other params...]
```

Now time to sit back and feel good for helping yourself and all institutional Nextflow users make their experience of running pipelines as smooth and efficient as possible 😎.

## Additional tips and tricks

### If your cluster often have too-full harddisk space issues

Some clusters often have very strict user HDD-space footprint restrictions (or just is often full).

You can minimise the footprint that Nextflow runs use by using the `cleanup` directive.

This directive goes _outside_ any of the dedicated scopes, and is simply defined as:

```nextflow
cleanup = true
```

This directive will, on a _successful_ completion of a Nextflow run, automatically delete all intermediate files stored in the `work/` directory. Note that none of the 'published' files in the `--outdir` results directory will be deleted if the pipeline specifies to _copy_ results files, _however_ removal of these intermediate files will mean debugging 'silent fails' more difficult.

To get around this, we suggest that you can optionally make an additional internal profile called debug:

```nextflow
debug {
  cleanup = false
}
```

which allows you to override the default `cleanup = true` behaviour, if you've set this as so.

Alternatively, if your cluster utilises a scratch space to store intermediate files, this can be specified a the `process` scope of the profile with the [`scratch`](https://www.nextflow.io/docs/latest/process.html#scratch) directive. For example:

```nextflow
process {
  <other_directives>
  scratch = '/<path>/<to>/<scratch>/<my_scratch>/'
}
```
