---
title: Troubleshooting
subtitle: How to troubleshoot common mistakes and issues
---

## Input files not found

If only no file, only one input file , or only read one and not read two is picked up then something is wrong with your input file declaration

1. The path must be enclosed in quotes (`'` or `"`)
2. The path must have at least one `*` wildcard character. This is even if you are only running one paired end sample.
3. When using the pipeline with paired end data, the path must use `{1,2}` or `{R1,R2}` notation to specify read pairs.
4. If you are running Single end data make sure to specify `--singleEnd`

If the pipeline can't find your files then you will get the following error

```bash
ERROR ~ Cannot find any reads matching: *{1,2}.fastq.gz
```

Note that if your sample name is "messy" then you have to be very particular with your glob specification. A file name like `L1-1-D-2h_S1_L002_R1_001.fastq.gz` can be difficult enough for a human to read. Specifying `*{1,2}*.gz` wont work give you what you want Whilst `*{R1,R2}*.gz` will.

## Data organization
The pipeline can't take a list of multiple input files - it takes a glob expression. If your input files are scattered in different paths then we recommend that you generate a directory with symlinked files. If running in paired end mode please make sure that your files are sensibly named so that they can be properly paired. See the previous point.

## Image cannot be build

Sometimes, `mksquashfs` cannot be found on the login node or workstation that you intend to use, thus the Singularity Image build fails unfortunately. See below code snippet that shows such a typical failure:

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

## Cannot find input files when using Singularity

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

## Extra resources and getting help

If you still have an issue with running the pipeline then feel free to contact us via the [Slack](https://nf-co.re/join/slack) channel or by opening an issue in the respective pipeline repository on GitHub asking for help.

If you have problems that are directly related to Nextflow and not our pipelines or the nf-core framework [tools](https://github.com/nf-core/tools) then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).
