---
title: Configuration
subtitle: How to troubleshoot common mistakes and issues
weight: 20
---

## Manual debugging on clusters using schedulers

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
