---
title: Troubleshooting
subtitle: How to troubleshoot common mistakes and issues
---

- [Troubleshooting basics](#troubleshooting-basics)
  - [Start small](#start-small)
  - [Categorize the type of error](#categorize-the-type-of-error)
  - [Read the log and check the work directory](#read-the-log-and-check-the-work-directory)
  - [Asking for help](#asking-for-help)
  - [Troubleshooting talk](#troubleshooting-talk)
- [Configuration](#configuration)
  - [Manual debugging on clusters using schedulers](#manual-debugging-on-clusters-using-schedulers)
- [Input files not found](#input-files-not-found)
  - [Direct input](#direct-input)
  - [Output for only a single sample although I specified multiple with wildcards](#output-for-only-a-single-sample-although-i-specified-multiple-with-wildcards)
  - [Sample sheet input](#sample-sheet-input)
  - [Data organization](#data-organization)
- [The pipeline crashes almost immediately with an early pipeline step](#the-pipeline-crashes-almost-immediately-with-an-early-pipeline-step)
  - [Tool not found](#tool-not-found)
  - [Error related to Docker](#error-related-to-docker)
  - [Error related to Singularity](#error-related-to-singularity)
  - [Error related to HPC Schedulers](#error-related-to-hpc-schedulers)
- [Cannot find input files when using Singularity](#cannot-find-input-files-when-using-singularity)
- [Warning about sticked on revision](#warning-about-sticked-on-revision)
- [My pipeline crashes part way through a run at a certain step with a non 0 exit code](#my-pipeline-crashes-part-way-through-a-run-at-a-certain-step-with-a-non-0-exit-code)
- [I get a exceeded job memory limit error](#i-get-a-exceeded-job-memory-limit-error)
- [Crashed pipeline with an error but Nextflow is still running](#crashed-pipeline-with-an-error-but-nextflow-is-still-running)
  - [A step of a pipeline wasn't executed](#a-step-of-a-pipeline-wasnt-executed)
- [My pipeline update doesn't seem to do anything](#my-pipeline-update-doesnt-seem-to-do-anything)
- [Unable to acquire lock error](#unable-to-acquire-lock-error)
- [Using a local version of iGenomes](#using-a-local-version-of-igenomes)
- [Extra resources and getting help](#extra-resources-and-getting-help)

## Troubleshooting basics

These are the recommended steps for troubleshooting a pipeline.

### Start small

Before using the pipeline with your data, make sure to run a test using:

```bash
nextflow run nf-core/<pipeline_name> -profile test,docker
```

If Docker is not installed, you can replace `docker` with `singularity` or `conda`, see
the [Getting Started](https://nf-co.re/docs/usage/introduction) tutorial for further information. If a test fails, it might indicate that
there is an issue with the installation or configuration of Nextflow or software management tool, rather than a pipeline error.

You might also want to check the following:

1. Nextflow is up to date. Use `nextflow self-update` to update a typical installation or `conda update nextflow` for a Bioconda installation.
2. There is enough disk space, this will avoid running out of space while you are running the pipeline.
3. Docker daemon is running (if you are using Docker to manage dependencies).

### Categorize the type of error

For this step you try to identify when the error occurs:

1. Before the first process: errors that occur before the first process might be related to an outdated version of Nextflow, updating to the newest version could help solving the issue. An example error is:

   ```bash
   N E X T F L O W  ~  version 0.27.3
   Launching `./main.nf` [prickly_snyder] - revision: bb0fa33a13
   ERROR ~ Unknown config attribute: projectDir -- check config file:
   nextflow.config
   null
   -- Check '.nextflow.log' file for details
   ```

2. During the first process: when an error appears during the first process it might indicate an issue with software dependencies, to specify how Nextflow should handle dependencies you need to select a [configuration profile](https://nf-co.re/docs/usage/configuration#basic-configuration-profiles). This type of error might also be related to a missing command required to run the pipeline. Example error:

   ```bash
   Command exit status:
     127
   Command output:
     (empty)
   Command error:
     .command.sh: line 3: rsem-prepare-reference: command not found
   Work dir:
     /home/lfaller/nextflow/rnaseq/work/f7/b6ef5a3f12f5efbf641f19046aca74
   Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`
   Unexpected error [AbortedException]
    -- Check script '/home/lfaller/.nextflow/assets/nf-core/rnaseq/./workflows/rnaseq.nf' at line: 603 or see '.nextflow.log' file for more details
   ```

3. During run: for errors that occur while a pipeline is running or generating outputs it might be helpful to check log files as explained below.

4. While generating outputs: if an expected process output is missing, Nextflow will fail with the message `Missing output file(s)`. Then the error message of that tool will be displayed. Example error:

   ```bash
   [6:16 PM] Error executing process > 'FASTQC (hct116_h3k4me1_IP_R1_T1)'
   Caused by:
     Missing output file(s) `*.{zip,html}` expected by process `FASTQC (hct116_h3k4me1_IP_R1_T1)`
   Command executed:
     [ ! -f  hct116_h3k4me1_IP_R1_T1.fastq.gz ] && ln -s hct116_h3k4me1_clean.fastq.gz
   hct116_h3k4me1_IP_R1_T1.fastq.gz
     fastqc -q -t 6 hct116_h3k4me1_IP_R1_T1.fastq.gz
   Command exit status:
     0
   Command output:
     (empty)
   Command error:
     WARNING: Your kernel does not support swap limit capabilities or the cgroup is not mounted.
   Memory limited without swap.
     Failed to process file hct116_h3k4me1_IP_R1_T1.fastq.gz
     uk.ac.babraham.FastQC.Sequence.SequenceFormatException: Ran out of data in the middle of a
   fastq entry.  Your file is probably truncated
      at uk.ac.babraham.FastQC.Sequence.FastQFile.readNext(FastQFile.java:179)
      at uk.ac.babraham.FastQC.Sequence.FastQFile.next(FastQFile.java:125)
      at uk.ac.babraham.FastQC.Analysis.AnalysisRunner.run(AnalysisRunner.java:77)
      at java.base/java.lang.Thread.run(Thread.java:834)
   ```

### Read the log and check the work directory

Checking the log files can help you to identify the type of error and where the error occurred. In order to search the output related to the error we need to understand the anatomy of the work directory:

1. `command.out` STDOUT from tool.
2. `command.err` STDERR from tool.
3. `command.log` contains both STDOUT and STDERR from tool.
4. `command.begin` created as soon as the job launches.
5. `exitcode` created when the job ends, with exit code.
6. `command.trace` logs of compute resource usage.
7. `command.run` wrapper script used to run the job.
8. `command.sh` process command used for this task.

If you checked the files and identified the type of error and where it occurred but were unable to solve it you can always ask for help.

### Asking for help

If you still have an issue with running the pipeline then feel free to contact us via the [Slack](https://nf-co.re/join/slack) channel. Please, consider the following guidelines:

- Pick the correct Slack channel to post in.
- Provide as much information as you can.
  - As a minimum the command and configs you used.
  - Use a Slack thread under your message if in doubt.
- Use markdown code blocks.
- Narrow the issue down as much as possible before asking.
- Explain the steps to reproduce if possible.

You can also open an issue in the respective pipeline repository on GitHub asking for help. In order to open the issue:

- Narrow the issue down as much as possible before opening the issue.
- Fill in the bug issue template.
- Explain the steps to reproduce.
- If you think you know the solution, please say so.
- If you think you can fix the problem, please make a pull request.

If you have problems that are directly related to Nextflow and not our pipelines or the nf-core framework [tools](https://github.com/nf-core/tools) then check out the [Nextflow Slack Channel](https://nextflow.io/slack-invite.html).

## Troubleshooting talk

For more detailed information about troubleshooting a pipeline run, you can check out the [Bytesize talk by Phil Ewels](https://www.youtube.com/embed/z9n2F4ByIkY) and the [accompanying slides](https://widgets.figshare.com/articles/19382933/embed?show_title=1).

## Configuration

### Manual debugging on clusters using schedulers

In some cases, when testing configuration files on clusters using a scheduler, you may get failed jobs with 'uninformative' errors and vague information being given for the cause.

In such cases, a good way of debugging such a failed job is to change to the working directory of the failed process (which should be reported by Nextflow), and try to _manually_ submit the job.

You can do this by submitting to your cluster the `.command.run` file found in the working directory using the relevant submission command.

For example, let's say you get an error like this on a SLURM cluster.

```console
Caused by:
  Failed to submit process to grid scheduler for execution

Command executed:

  sbatch .command.run

Command exit status:
  -

Command output:
  (empty)

Work dir:
  /<path>/<to>/work/e5/6cc8991c2b16c11a6356028228377e
```

This does not tell you _why_ the job failed to submit, but is often is due to a 'invalid' resource submission request, and the scheduler blocks it. But unfortunately, Nextflow does not pick the message reported by the cluster.

Therefore, in this case I would switch to the working directory, and submit the `.command.run` file using SLURM's `sbatch` command (for submitting batch scripts).

```bash
$ cd  /<path>/<to>/work/e5/6cc8991c2b16c11a6356028228377e
$ sbatch .command.run
sbatch: error: job memory limit for shared nodes exceeded. Must be <= 120000 MB
sbatch: error: Batch job submission failed: Invalid feature specification
```

In this case, SLURM has printed to my console the reason _why_, the job failed to be submitted.

With this information, I can go back to my configuration file, and tweak the settings accordingly, and run the pipeline again.

## Input files not found

If the pipeline can't find your files then you will get the following error

```bash
ERROR ~ Cannot find any reads matching: *{1,2}.fastq.gz
```

### Direct input

Or when you're using a input method like `--input '/<path>/<to>/*_fq.gz'`, but only pick up one file, or only one file per pair being processed during the run, please note the following:

1. [The path must be enclosed in quotes (`'` or `"`)](#output-for-only-a-single-sample-although-i-specified-multiple-with-wildcards)
2. The path must have at least one `*` wildcard character i.e. following a ['glob' pattern](<https://en.wikipedia.org/wiki/Glob_(programming)>). This is even if you are only running one paired end sample.
   - A description of valid pattern matching can be seen [here](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) for java and [here](https://www.nextflow.io/docs/latest/channel.html?highlight=glob#frompath) for Nextflow
3. When using the pipeline with paired end data, the path must use `{1,2}` or `{R1,R2}` notation to specify read pairs.
   - This notation is interpreted by Nextflow to mean anything with the same string other than R1 and R2 in the file name, will be be assumed to be a pair of files.
4. If you are running single-end data make sure to specify `--singleEnd`
5. [Your data should be organised in a 'tidy' manner](#data-organization)

A few examples are as follows:

- Running with a single, single-end FASTQ file as input (this will produce output files for this sample only)

  ```bash
  nextflow run nf-core/<pipeline> -input 'my_data.fastq.gz` --single_end
  ```

- Running multiple single-end FASTQ files as input using a wildcard glob pattern. This will find all files in the directory beginning with `my_`, and ending in `.fastq.gz`, with each file with any other characters between those two strings being considered distinct samples (and will produce output files for each of the multiple input files).

  ```bash
  nextflow run nf-core/<pipeline> -input 'my_*.fastq.gz` --single_end
  ```

- Running multiple paired-end FASTQ files as input using wildcard and grouping glob patterns. This will find all files in the directory beginning with `my_`, and ending in `.fastq.gz`, with each file with any other characters between those two strings being considered distinct samples. However, any pair of file names that are exactly the same other than `R1` and `R2` will be grouped together, and processed as related files. You will in most cases get output files for each distinct file, but with the `*{R1,R2}` syntax, R1 and R2 pairs are collapsed into one.

  ```bash
  nextflow run nf-core/<pipeline> -input 'my_*{R1,R2}.fastq.gz`
  ```

Note that if your sample name is "messy" then you have to be very particular with your glob specification (see point 2 above). A file name like `L1-1-D-2h_S1_L002_R1_001.fastq.gz` can be difficult enough for a human to read. Specifying `*{1,2}*.gz` will not give you what you want, whilst `*{R1,R2}*.gz` will.

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

### Sample sheet input

If you are using a sample sheet or TSV input method, check there is not a mistake or typo in the path in a given column. Common mistakes are a trailing space at the end of the path, which can cause problems.

### Data organization

The pipeline can't take a list of multiple input files - it takes a glob expression. If your input files are scattered in different paths then we recommend that you generate a directory with symlinked files. If running in paired end mode please make sure that your files are sensibly named so that they can be properly paired. See the previous point.

## The pipeline crashes almost immediately with an early pipeline step

Sometimes a newly downloaded and set up nf-core pipeline will encounter an issue where a run almost immediately crashes (e.g. at `fastqc`, `output_documentation` etc.) saying the tool could not be found or similar.

The first thing to do is always check the `.nextflow.log` to see if it reports contains specific error. Common cases are described below.

:warning: Note that just because Nextflow reports a particular tool failed, this _does not_ necessarily mean it's an issue with the tool itself. It's important to always _fully_ read the error message to identify possible causes.

### Tool not found

The most common case is when a user has forgotten to specify a container/environment profile.

If you do not specify one with the `-profile`, Nextflow by default looks for all the required tools that will be all manually installed on your machine/cluster and specified `$PATH`.

It is _not_ recommended to run without a container/environment system as then your analysis will not be reproducible by others. You should pick one of: `docker`, `singularity`, `podman` and `conda` - depending on which your system already has available. See the nf-core [Installation](https://nf-co.re/docs/usage/installation) documentation for more information.

### Error related to Docker

You may have an outdated container. This happens more often when running on the `dev` branch of a nf-core pipeline, because Docker will _not_ update the container on each new commit, and thus may not get new tools called within the pipeline code.

To fix, just re-pull the pipeline's Docker container manually with:

```bash
docker pull nfcore/<pipeline>:dev
```

### Error related to Singularity

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

### Error related to HPC Schedulers

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

This is a Nextflow error, with less-commonly seen Git 'terminology'. What this means is that you have multiple versions of the pipeline pulled (e.g. 2.0.0, 2.1.0, 2.1.1, dev etc.), and it is not sure which one to use. Therefore, with every `nextflow run nf-core/<PIPELINE>` command you should always indicate which version with `-r`.

For example:

```bash
nextflow run nf-core/<pipeline> -r 2.1.0 --input '/<path>/<to>/data/*_{R1,R2}_*.fq.gz' <...>
```

Specifying the version of the run you are using is highly recommended, as it helps in full reproducibility. In the sense that if you explicitly record the whole command _with_ the version for your publication or internal reports, then anyone who wants to check your work can use the exact version you used (including all internal tools).

You can see more information on the Nextflow documentation [here](https://www.nextflow.io/docs/latest/sharing.html#handling-revisions).

## My pipeline crashes part way through a run at a certain step with a non 0 exit code

Sometimes part way through a run, a particular tool or step of the pipeline will fail. While Nextflow and nf-core pipelines try to solve some of these issues for you, this is not always possible. If the particular eventually tool fails, Nextflow will report the process that failed with a non-0 error code, and print the command that failed.

An example is as follows:

```text
Error executing process > 'markduplicates (DA117)'

Caused by:
  Process `markduplicates (DA117)` terminated with an error exit status (137)

Command executed:

  picard -Xmx16384M -Xms16384M MarkDuplicates INPUT=DA117.bam OUTPUT=DA117_rmdup.bam REMOVE_DUPLICATES=TRUE AS=TRUE METRICS_FILE="DA117_rmdup.metrics" VALIDATION_STRINGENCY=SILENT
  samtools index DA117_rmdup.bam

Command exit status:
  137
```

> :warning: Each exit code can mean different things to different tools as well as in different environments. Therefore it is not always easy for developers to predict the exact issue and solution!

Common exit codes and and **_potential_** solutions are as follows:

| Exit Code | Possible Cause | Solution                                                                                                                                                                                                     |
| --------- | -------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `104`     | out of memory  | increase memory of process or number of retries in profile: [Quick reference](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), [Step By Step](#i-get-a-exceeded-job-memory-limit-error) |
| `134`     | out of memory  | increase memory of process or number of retries in profile: [Quick reference](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), [Step By Step](#i-get-a-exceeded-job-memory-limit-error) |
| `137`     | out of memory  | increase memory of process or number of retries in profile: [Quick reference](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), [Step By Step](#i-get-a-exceeded-job-memory-limit-error) |
| `139`     | out of memory  | increase memory of process or number of retries in profile: [Quick reference](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), [Step By Step](#i-get-a-exceeded-job-memory-limit-error) |
| `143`     | out of memory  | increase memory of process or number of retries in profile: [Quick reference](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), [Step By Step](#i-get-a-exceeded-job-memory-limit-error) |
| `247`     | out of memory  | increase memory of process or number of retries in profile: [Quick reference](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), [Step By Step](#i-get-a-exceeded-job-memory-limit-error) |

If in doubt, Google is your friend! Many exit codes are roughly similar across many tools; even if search results don't mention your tool exactly, you can try a solution similar to the one proposed for the other tool.

If you are still unable to resolve the issue, please make a GitHub issue on the corresponding pipeline repository.

## I get a exceeded job memory limit error

While Nextflow tries to make your life easier by automatically retrying jobs that run out of memory with more resources (until a [specified max-limit](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources)), sometimes you may have such large data that you run out even after the default 3 retries.

To fix this you need to change the default memory requirements for the process that is breaking. We can do this by making a custom profile, which we then provide to the Nextflow run command.

For example, let's say it's the `markduplicates` process that is running out of memory (as displayed on the Nextflow running display).

- First we need to check to see what default memory value we have. Go to the main code of the pipeline by going to the corresponding pipeline GitHub repository and open the `main.nf` file. Search in your browser for `process markduplicates`.
- Once found, check the line called `label` and note down the corresponding label. In this case the label could be `process_low`.
- Go back to the main github repository, and open `conf/base.config`. Now, search in your browser for `withLabel:'process_low'`.
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

  - Note that with the above example you will **_not_** have the automatic retry mechanism that resubmits jobs with increased resource requests (given appropriate exit codes). The job will still be resubmitted on failure but with `16.GB` each time.

    - If you want this, use the following syntax instead:

      ```nextflow
      memory = { check_max( 16.GB * task.attempt, 'memory' ) }
      ```

    - Next, copy the `check_max()` function from the pipeline's `nextflow.config` file (e.g. [here](https://github.com/nf-core/rnaseq/blob/3643a94411b65f42bce5357c5015603099556ad9/nextflow.config#L190-L221)) to the bottom of your custom config file.
    - `16.GB * task.attempt` multiplies the memory request by the index of the retry. So if the job failed and is being tried a second time, it requests `32.GB`.
    - The `check_max()` function prevents Nextflow requesting excessive resources above what is available on your system. This effectively sets a ceiling on the resources and prevents the pipeline from crashing if it goes too high. Unfortunately because of the order in which pipeline code and Nextflow configs are parsed, this function needs to be defined in your custom config file.

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

## Using a local version of iGenomes

The iGenomes config file uses `params.igenomes_base` to make it possible to use a local copy of iGenomes. However, as custom config files are loaded after `nextflow.conf` and the `igenomes.config` has already been imported and parsed, setting `params.igenomes_base` in a custom config file has no effect and the pipeline will use the s3 locations by default.

You can specify a local iGenomes path by either:

- Setting the igenomes_base path in a configuration profile.

```nextflow
params {
    igenomes_base = '/<path>/<to>/<data>/igenomes'
}
```

- Specifying an `--igenomes_base` path in your execution command.

```bash
nextflow run nf-core/<pipeline> --input <input> -c <config> -profile <profile> --igenoms_base <path>/<to>/<data>/igenomes
```

- Specifying the `igenomes_base` parameter in a `params` file provided with `-params-file` in `yaml` or `json` format.

```bash
nextflow run nf-core/<pipeline> -profile <profile> -params-file params.yml
```

Where the `params.yml` file contains the pipeline params:

```bash
input: '/<path>/<to>/<data>/input'
igenomes_base: '/<path>/<to>/<data>/igenomes'
```

## Extra resources and getting help

If you still have an issue with running the pipeline then feel free to contact us via the [Slack](https://nf-co.re/join/slack) channel or by opening an issue in the respective pipeline repository on GitHub asking for help.

If you have problems that are directly related to Nextflow and not our pipelines or the nf-core framework [tools](https://github.com/nf-core/tools) then checkout the [Slack community chat](https://nextflow.io/slack-invite.html) or the [discussion forum](https://github.com/nextflow-io/nextflow/discussions) on GitHub.
