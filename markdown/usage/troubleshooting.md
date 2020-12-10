---
title: Troubleshooting
subtitle: How to troubleshoot common mistakes and issues
---

<!-- TOC -->

- [Input files not found](#input-files-not-found)
  - [Direct input](#direct-input)
  - [Sample sheet input](#sample-sheet-input)
  - [Output for only a single sample although I specified multiple with wildcards](#output-for-only-a-single-sample-although-i-specified-multiple-with-wildcards)
  - [Data organization](#data-organization)
- [The pipeline crashes almost immediately with an early pipeline step](#the-pipeline-crashes-almost-immediately-with-an-early-pipeline-step)
  - [Error related to Docker](#error-related-to-docker)
  - [Error related to Singularity](#error-related-to-singularity)
- [Cannot find input files when using Singularity](#cannot-find-input-files-when-using-singularity)
- [Warning about sticked on revision](#warning-about-sticked-on-revision)
- [I get a exceeded job memory limit error](#i-get-a-exceeded-job-memory-limit-error)
- [Crashed pipeline with an error but Nextflow is still running](#crashed-pipeline-with-an-error-but-nextflow-is-still-running)
  - [A step of a pipeline wasn't executed](#a-step-of-a-pipeline-wasnt-executed)
- [My pipeline update doesn't seem to do anything](#my-pipeline-update-doesnt-seem-to-do-anything)
- [Unable to acquire lock error](#unable-to-acquire-lock-error)
- [Extra resources and getting help](#extra-resources-and-getting-help)

<!-- /TOC -->

## Input files not found

If the pipeline can't find your files then you will get the following error

```bash
ERROR ~ Cannot find any reads matching: *{1,2}.fastq.gz
```

### Direct input

Or when you're using a input method like `--input '/<path>/<to>/*_fq.gz'`, but only pick up one file, or only one file per pair being processed during the run, please note the following:

1. [The path must be enclosed in quotes (`'` or `"`)](#output-for-only-a-single-sample-although-i-specified-multiple-with-wildcards)
2. The path must have at least one `*` wildcard character. This is even if you are only running one paired end sample.
3. When using the pipeline with paired end data, the path must use `{1,2}` or `{R1,R2}` notation to specify read pairs.
    - This notation is interpreted by Nextflow to mean anything with the same string other than R1 and R2 in the file name, will be be assumed to be a pair of files.
4. If you are running single-end data make sure to specify `--singleEnd`
5. [Your data should be organised in a 'tidy' manner](#data-organization)

Note that if your sample name is "messy" then you have to be very particular with your glob specification. A file name like `L1-1-D-2h_S1_L002_R1_001.fastq.gz` can be difficult enough for a human to read. Specifying `*{1,2}*.gz` wont work give you what you want, whilst `*{R1,R2}*.gz` will.

### Sample sheet input

If you are using a sample sheet or TSV input method, check there is not a mistake or typo in the path in a given column. Common mistakes ar e a trailing space at the end of the path, which can cause problems.

### Output for only a single sample although I specified multiple with wildcards

You must specify paths to files in quotes, otherwise your _shell_ (e.g. bash) will evaluate any wildcards (\*) rather than Nextflow.

For example:

```bash
nextflow run nf-core/<pipeline> --input /path/to/sample_*/*.fq.gz
```

Maybe evaluated by your shell as:

```bash
nextflow run nf-core/<pipeline> --input /path/to/sample_1/sample_1.fq.gz /path/to/sample_1/sample_1.fq.gz /path/to/sample_1/sample_1.fq.gz
```

And Nextflow will only take the first path after `--input`, ignoring the others.

On the other hand, encapsulating the path in quotes will allow _Nextflow_ to evaluate the paths.

```bash
nextflow run nf-core/<pipeline> --input "/path/to/sample_*/*.fq.gz"
```

### Data organization

The pipeline can't take a list of multiple input files - it takes a glob expression. If your input files are scattered in different paths then we recommend that you generate a directory with symlinked files. If running in paired end mode please make sure that your files are sensibly named so that they can be properly paired. See the previous point.

## The pipeline crashes almost immediately with an early pipeline step

Sometimes a newly downloaded and set up nf-core pipeline will encounter an issue where a run almost immediately crashes (e.g. at `fastqc`, `output_documentation` etc.) saying the tool could not be found or similar.

The first thing to do is always check the `.nextflow.log` to see if it reports contains specific error. Common cases are described below.

### Error related to Docker

You may have an outdated container. This happens more often when running on the `dev` branch of a nf-core pipeline, because Docker will _not_ update the container on each new commit, and thus may not get new tools called within the pipeline code.

To fix, just re-pull the pipeline's Docker container manually with:

```bash
docker pull nfcore/<pipeline>:dev
```

### Error related to Singularity

If you're running Singularity, it could be that Nextflow cannot access your Singularity image properly - often due to missing bind paths. See [_Cannot find input files when using Singularity_](https://nf-co.re/usage/troubleshooting#cannot-find-input-files-when-using-singularity) for more information.

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

## Warning about sticked on revision

If you get a warning like the following:

```bash
Project nf-core/<pipeline> currently is sticked on revision: dev -- you need to specify explicitly a revision with the option -r to use it
```

This is a Nextflow error, with less commonly seen terminology. What this means is that you have multiple versions of the pipeline pulled (e.g. 2.0.0, 2.1.0, 2.1.1, dev etc.). When you have multiple versions, you must always specify which one you want to use - there is no 'default'. Therefore, with every `nextflow run nf-core/<PIPELINE>` command you must always indicate which version with `-r`.

For example:

```bash
nextflow run nf-core/<pipeline> -r 2.1.0 --input '/<path>/<to>/data/*_{R1,R2}_*.fq.gz' <...>
```

Specifying the version of the run you are using is highly recommended, as it helps in full reproducibility. In the sense that if you explicitly record the whole command _with_ the version for your publication or internal reports, then anyone who wants to check your work can use the exact version you used (including all internal tools).

You can see more information on the Nextflow documentation [here](https://www.nextflow.io/docs/latest/sharing.html#handling-revisions).

## I get a exceeded job memory limit error

While Nextflow tries to make your life easier by automatically retrying jobs that run out of memory with more resources (until a [specified max-limit](https://nf-co.re/usage/configuration#tuning-workflow-resources)), sometimes you may have such large data that you run out even after the default 3 retries.

To fix this you need to change the default memory requirements for the process that is breaking. We can do this by making a custom profile, which we then provide to the Nextflow run command.

For example, lets say it's the `markduplicates` process that is running out of memory (as displayed on the Nextflow running display).

- First we need to check to see what default memory value we have. Go to the main code of the pipeline by going to the corresponding pipeline GitHub repository and open the `main.nf` file. Use your browser's find functionality for: `process markduplicates`.
- Once found, check the line called `label` and note down the corresponding label. In this case the label could be `process_low`.
- Go back to the main github repository, and open `conf/base.config`. Again use the browser;s find functionality to search for: `withLabel:'process_low'`.
- Note what the `memory` field is set to (e.g. `4.GB`) on a line like: `memory = { check_max( 4.GB * task.attempt, 'memory' )})`.
- Back on your working machine, make a new text file called `custom_resources.conf`. This should be saved somewhere centrally so you can reuse it.
    > If you think this would be useful for multiple people in your lab/institute, we highly recommend you make an institutional profile at [nf-core/configs](https://github.com/nf-core/configs). This will simplify this process in the future.
- Within this file, add the following. Note we have increased the default `4.GB` to `16.GB`.

    ```nextflow
    profiles {
        big_data {
          process {
            withName: markduplicates {
              memory = 16.GB
            }
          }
        }
    }
    ```

    > Note that with this you will _not_ have the automatic retry mechanism. If you want this, re-add the `check_max()` function on the `memory` line (as in the pipeline's `conf/base.config`) file, and add to the bottom of the entire file (outside the profiles block!). This block is the one starting with `def check_max(obj, type) {`, which is at the end of the pipeline's `nextflow.config` file.

- Once saved, modify the original Nextflow run command:

    ```bash
    nextflow run nf-core/<pipeline> -c /<path>/<to>/custom_resources.conf -profile big_data,<original>,<profiles> <...>
    ```

  - We have added `-c` to specify which file to use for the custom profiles, and then added the `big_data` profile to the original profiles you were using.
  - :warning: it's important that the `big_data` profile name comes first, to ensure it overwrites any parameters set in the subsequent profiles. Profile names should be comma separated with no spaces.

## Crashed pipeline with an error but Nextflow is still running

If this happens, you can either wait until all other already running jobs to safely finish, or if Nextflow _still_ does not stop press `ctrl + c` on your keyboard (or equivalent) to stop the Nextflow run.

> :warning: if you do this, and do not plan to fix the run make sure to delete the `work` folder generated that is generated at the same as `results` (or specified with the Nextflow variable `-w`). Otherwise you may end up a lot of large intermediate files being left! You can clean a Nextflow run of all intermediate files with `nextflow clean -f -k` or delete the `work/` directory.

### A step of a pipeline wasn't executed

Possible options:

1. If an optional step, check for a typo in the parameter name. Nextflow _does not_ check for this
2. Check that an upstream step/process was turned on (if a step/process requires the output of an earlier process, it will not be activated unless it receives the output of that process)

## My pipeline update doesn't seem to do anything

To download a new version of a pipeline, you can use the following, replacing `<version>` to the corresponding version.

```bash
nextflow pull nf-core/<pipeline> -r <version>
```

However, in very rare cases, minor fixes to a version will be pushed out without a version number bump. This can confuse Nextflow slightly, as it thinks you already have the 'broken' version from your original pipeline download.

> This _shouldn't_ happen with stable versions and normally only happens on `dev` branches.

If when running the pipeline you don't see any changes in the fixed version when running it, you can try removing your Nextflow's nf-core pipeline cache typically stored in your home directory with:

```bash
rm -r ~/.nextflow/assets/nf-core/<pipeline>
```

And re-pull the pipeline with the command above. This will install a fresh version of the version with the fixes.

## Unable to acquire lock error

Errors like the following:

```bash
Unable to acquire lock on session with ID 84333844-66e3-4846-a664-b446d070f775
```

Normally suggest a previous Nextflow run (on the same folder) was not cleanly killed by a user (e.g. using ctrl + z to hard kill a crashed run).

To fix this, you must clean the entirety of the run's `work/` directory e.g. with `rm -r work/` and re-running from scratch.

`ctrl +z` is **not** a recommended way of killing a Nextflow job. Runs that take a long time to fail are often still running because other job submissions are still running. Nextflow will normally wait for those processes to complete before cleaning shutting down the run (to allow rerunning of a run with `-resume`). `ctrl + c` is much safer as it will tell Nextflow to stop earlier but cleanly.

## Extra resources and getting help

If you still have an issue with running the pipeline then feel free to contact us via the [Slack](https://nf-co.re/join/slack) channel or by opening an issue in the respective pipeline repository on GitHub asking for help.

If you have problems that are directly related to Nextflow and not our pipelines or the nf-core framework [tools](https://github.com/nf-core/tools) then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).
