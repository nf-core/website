---
title: Early crash
subtitle: How to troubleshoot common mistakes and issues
---

### The pipeline crashes almost immediately with an early pipeline step

Sometimes a newly downloaded and set up nf-core pipeline will encounter an issue where a run almost immediately crashes (e.g. at `fastqc`, `output_documentation` etc.) saying the tool could not be found or similar.

The first thing to do is always check the `.nextflow.log` to see if it reports contains specific error. Common cases are described below.

:::warning
Note that just because Nextflow reports a particular tool failed, this _does not_ necessarily mean it's an issue with the tool itself. It's important to always _fully_ read the error message to identify possible causes.
:::

#### Tool not found

The most common case is when a user has forgotten to specify a container/environment profile.

If you do not specify one with the `-profile`, Nextflow by default looks for all the required tools that will be all manually installed on your machine/cluster and specified `$PATH`.

It is _not_ recommended to run without a container/environment system as then your analysis will not be reproducible by others. You should pick one of: `docker`, `singularity`, `podman` and `conda` - depending on which your system already has available. See the nf-core [Installation](https://nf-co.re/docs/usage/installation/) documentation for more information.

#### Error related to Docker

You may have an outdated container. This happens more often when running on the `dev` branch of a nf-core pipeline, because Docker will _not_ update the container on each new commit, and thus may not get new tools called within the pipeline code.

To fix, just re-pull the pipeline's Docker container manually with:

```bash
docker pull nfcore/<pipeline>:dev
```

If you work in a `dev` branch, you may also want to consider putting a request for a pull in each run you do with nextflow, by putting this line of code in your nextflow.config file:

```groovy
docker {
    enabled = true
    runOptions = '--pull=always'
}
```

#### Error related to Singularity

If you're running Singularity, it could be that Nextflow cannot access your Singularity image properly - often due to missing bind paths. See [_Cannot find input files when using Singularity_](https://nf-co.re/docs/usage/troubleshooting#cannot-find-input-files-when-using-singularity) for more information.

Sometimes, `mksquashfs` cannot be found on the login node or workstation that you intend to use, thus the Singularity image build fails unfortunately. See below code snippet that shows such a typical failure:

```bash
Caused by:
  Failed to pull singularity image
  command: singularity pull --name nfcore-rnaseq-1.3.img docker://nfcore/rnaseq:1.3 > /dev/null
  status : 255
  .....
    INFO:    Creating SIF file...
    FATAL:   Unable to pull docker://nfcore/rnaseq:1.3: While searching for mksquashfs: exec: "mksquashfs": executable file not found in $PATH
```

If this is the case, please install `mksquashfs` or ask your IT department to install the package for you.

#### Error related to HPC Schedulers

If working on a cluster, pipelines can crash if the profile used is not correctly configured for that environment. Typical issues can include missing cluster profile in `-profile`, incorrectly specified executor, or incompatible memory/CPU node maximums set in that institutional profile. See [nf-core/configs](https://github.com/nf-core/configs) and [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#) for more information.

These types of error can look like the following:

```bash
Caused by:
  Failed to submit process to grid scheduler for execution

Command executed:

  sbatch .command.run

Command exit status:
  1

Command output:
  sbatch: error: Batch job submission failed: Invalid account or account/partition combination specified
```

#### Cannot find input files when using Singularity

Depending on how you install Singularity on your system several reoccurring issues have been reported. Typically these result in error messages such as these:

```bash

Command error:
  �[91mERROR  : Failed to resolve path to /home/bla/singularity/mnt/container: No such file or directory
  �[0m�[31mABORT  : Retval = 255
  �[0m
```

You can verify if this is indeed a problem by using a Singularity Shell to access your container, and to check whether the required paths are available **inside** your container:

```bash
singularity shell /path/to/work/singularity/imagename.simg
ls -la /path/to/work
ls -la /path/to/network_storage
```

If any of these `ls -la` commands returns a `Not Found` error, you will need to set/update Singularity Bind Paths on your system.

The Singularity installation requires certain host paths to be bound. Please see [Singularity Bind Paths Documentation](https://sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html) for a more detailed explanation. In many cases this can be resolved by adding these paths to your `/etc/singularity/singularity.conf` as highlighted in the documentation:

```bash
bind path = /beegfs/work/
bind path = /scratch
bind path = /gpfs
bind path = /home
```

Alternatively, you can also add Singularity Bind Paths to your Nextflow call, e.g. using `autoMounts` and/or `runOptions` in the [Singularity scope](https://www.nextflow.io/docs/latest/config.html#config-singularity)
