---
title: Crash halfway through
subtitle: How to troubleshoot common mistakes and issues
---

### Crashes at a certain step with a non 0 exit code

Sometimes part way through a run, a particular tool or step of the pipeline will fail. While Nextflow and nf-core pipelines try to solve some of these issues for you, this is not always possible. If the particular eventually tool fails, Nextflow will report the process that failed with a non-0 error code, and print the command that failed.

An example is as follows:

```text
Error executing process > 'NFCORE_SAREK:SAREK:MARKDUPLICATES (DA117)'

Caused by:
  Process `NFCORE_SAREK:SAREK:MARKDUPLICATES (DA117)` terminated with an error exit status (137)

Command executed:

  picard -Xmx16384M -Xms16384M MarkDuplicates INPUT=DA117.bam OUTPUT=DA117_rmdup.bam REMOVE_DUPLICATES=TRUE AS=TRUE METRICS_FILE="DA117_rmdup.metrics" VALIDATION_STRINGENCY=SILENT
  samtools index DA117_rmdup.bam

Command exit status:
  137
```

:::warning
Each exit code can mean different things to different tools as well as in different environments. Therefore it is not always easy for developers to predict the exact issue and solution!
:::

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

### I get a exceeded job memory limit error

If you hit an error such as `Process requirement exceeds available memory -- req: 100 GB; avail: 15.5 GB`, this implies the pipeline run has requested more memory than is available on your system.

While Nextflow tries to make your life easier by automatically retrying jobs that run out of memory with more resources (until a [specified max-limit](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources)), sometimes you may have such large data that you run out even after the default 3 retries.

To address this you will need to make a configuration file that tells Nextflow how much memory or CPUs are available on your system, and also how much memory or CPUs the failing step should use.

Please see the two following sections of the [configuring nf-core pipelines]
(https://nf-co.re/docs/usage/configuration) page:

- Setting limits: [Max Resources](https://nf-co.re/docs/usage/configuration#max-resources)
- Adusting the process resources: [Tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources)

The resulting configuration file can then be passed to your Nextflow run command with `-c <config_file> -resume` to resume the failed run but with the updated resource requirements.

### Crashed pipeline with an error but Nextflow is still running

If this happens, you can either wait until all other already running jobs to safely finish, or if Nextflow _still_ does not stop press `ctrl + c` on your keyboard (or equivalent) to stop the Nextflow run.

:::warning
if you do this, and do not plan to fix the run make sure to delete the `work` folder generated that is generated at the same as `results` (or specified with the Nextflow variable `-w`). Otherwise you may end up a lot of large intermediate files being left! You can clean a Nextflow run of all intermediate files with `nextflow clean -f -k` or delete the `work/` directory.
:::

#### A step of a pipeline wasn't executed

Possible options:

1. If an optional step, check for a typo in the parameter name. Nextflow _does not_ check for this (unless you created a workflow with the nf-core template, which provides parameter validation)
2. Check that an upstream step/process was turned on (if a step/process requires the output of an earlier process, it will not be activated unless it receives the output of that process)
