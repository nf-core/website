---
title: Troubleshooting
subtitle: Troubleshoot nf-core pipeline runs
shortTitle: Troubleshooting
---

This page covers common issues you may encounter when running nf-core pipelines.

## Running pipelines

### Pipeline crashes immediately

**Problem:**

Pipeline fails at the first process.

**Cause:**

No container or environment profile specified.

**Solution:**

Add `-profile docker`, `-profile singularity`, or `-profile conda` to your command:

```bash
nextflow run nf-core/pipeline -profile docker
```

### Cluster job submission fails

**Error message:**

```console
Failed to submit process to grid scheduler for execution.
```

**Solution:**

Manually submit the failed job to see the actual error:

1. Navigate to the work directory shown in the error
2. Submit the command file directly to your scheduler:

   ```bash
   sbatch .command.run  # For Slurm
   qsub .command.run    # For SGE
   ```

3. The scheduler will show the real error (for example, "job memory limit exceeded")
4. Adjust your configuration based on the error message

### Invalid executor configuration

**Error message:**

```console
sbatch: error: Invalid account or account/partition combination specified.
```

**Cause:**

- Missing cluster profile in `-profile` parameter
- Incorrectly specified executor in your configuration
- Resource requests exceed your cluster limits

**Solution:**

- Check if your institution has a profile in [nf-core/configs](https://github.com/nf-core/configs)
- Verify your executor name matches your scheduler (for example, `slurm` not `SLURM`)
- Confirm your `resourceLimits` match your cluster's capabilities

### Singularity bind path errors

**Error message:**

```console
ERROR: Failed to resolve path to /home/path: No such file or directory
```

**Cause:**

Singularity cannot access your file system paths.

**Solution:**

Add bind paths to your Nextflow configuration:

```groovy
singularity {
  enabled = true
  autoMounts = true
  runOptions = '-B /scratch -B /gpfs -B /work'
}
```

Or update your Singularity system configuration at `/etc/singularity/singularity.conf`.

### Container not updating

**Problem:**

Pipeline uses old tool versions despite specifying `dev` branch.

**Cause:**

Docker does not automatically update local images.

**Solution:**

Manually pull the latest container:

```bash
docker pull nfcore/pipeline:dev
```

Or add to your configuration to always pull:

```groovy
docker {
  enabled = true
  runOptions = '--pull=always'
}
```
