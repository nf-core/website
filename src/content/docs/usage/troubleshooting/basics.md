---
title: Troubleshooting basics
subtitle: How to troubleshoot common mistakes and issues
weight: 10
---

These are the recommended steps for troubleshooting a pipeline.

## Check Nextflow is definitely working

Before using the pipeline with your own data and parameters, make sure to run a test in a separate directory using:

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

## Categorize the type of error

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

## Read the log and check the work directory

Checking the log files can help you to identify the type of error and where the error occurred.

The first log file to look at is `.nextflow.log` that is placed in you executed `nextflow run`.
This reports all logging information from the overarching pipeline run.
Note this can be overwhelming! If so, proceed to the next step.

In order to search the output related to the error we need to understand the anatomy of the work directory. The work directory is reported at the end of the error message:

```console
Work dir:
  /home/harshil/repos/nf-core/fetchngs/work/c3/29eddb3ef99977c62f462b61f09afa
```

And contains:

1. `command.out` STDOUT from tool.
2. `command.err` STDERR from tool.
3. `command.log` contains both STDOUT and STDERR from tool.
4. `command.begin` created as soon as the job launches.
5. `exitcode` created when the job ends, with exit code.
6. `command.trace` logs of compute resource usage.
7. `command.run` wrapper script used to run the job.
8. `command.sh` process command used for this task.

If you checked the files and identified the type of error and where it occurred but were unable to solve it you can always ask for help.

## Asking for help

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
