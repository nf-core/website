---
title: Testing profiles
subtitle: Learn how to test your institutional profile before submission
shortTitle: Testing
weight: 6
---

Before submitting your institutional profile to nf-core/configs, test it thoroughly to verify that all settings work correctly. Testing prevents issues for future users and ensures smooth integration with nf-core pipelines.

## Fork and clone nf-core/configs

Test your profile using your fork of the nf-core/configs repository:

1. Fork the [nf-core/configs](https://github.com/nf-core/configs) repository to your GitHub account
1. Clone your fork to your local system or cluster:

   ```bash
   git clone https://github.com/<your_username>/configs.git
   cd configs
   ```

1. Create a new branch for your profile:

   ```bash
   git checkout -b add-big-university-profile
   ```

1. Add your configuration and documentation files as described in the file structure guide

## Test with custom config base

Use the `--custom_config_base` parameter to test your profile before it merges into the main repository. This parameter tells Nextflow to load configurations from your fork instead of the official nf-core/configs repository.

### Test from local files

If you cloned your fork to the cluster where you will run tests:

```bash
nextflow run nf-core/<pipeline_name> \
  -profile <your_profile>,test \
  --custom_config_base /path/to/your/configs
```

Replace:

- `<pipeline_name>`: The pipeline to test (for example, `rnaseq`, `sarek`, `atacseq`)
- `<your_profile>`: Your profile name (for example, `big_university`)
- `/path/to/your/configs`: The absolute path to your cloned configs directory

### Test from GitHub

If your changes are pushed to GitHub, test directly from your fork:

```bash
nextflow run nf-core/<pipeline_name> \
  -profile <your_profile>,test \
  --custom_config_base 'https://raw.githubusercontent.com/<your_username>/configs/<branch>'
```

Replace:

- `<your_username>`: Your GitHub username
- `<branch>`: Your branch name (for example, `add-big-university-profile`)

This approach lets you test without cloning to the cluster, but requires pushing changes to GitHub first.

## Use test profiles

Always test with a pipeline's built-in test profile to minimize runtime and resource usage. Test profiles use small datasets that complete quickly while exercising all pipeline features.

```bash
nextflow run nf-core/rnaseq \
  -profile big_university,test \
  --custom_config_base /path/to/your/configs
```

The `test` profile comes after your institutional profile so test settings override any conflicting parameters.

:::tip
Start with simple pipelines like `rnaseq` or `fetchngs` for initial testing. These pipelines have well-maintained test profiles and run relatively quickly.
:::

## What to verify

Check these aspects during testing:

### Profile loading

Verify your profile loaded correctly by checking the pipeline startup output:

```
-[nf-core/rnaseq] Pipeline Release: 3.14.0
-
-[nf-core/rnaseq] Config Profile Description: Big University HPC cluster profile
-[nf-core/rnaseq] Config Profile Contact: Jane Smith (@jsmith)
-[nf-core/rnaseq] Config Profile URL: https://hpc.biguniversity.edu
```

If you do not see your `config_profile_*` parameters, your profile did not load correctly. Check for syntax errors in your configuration file.

### Container engine detection

Confirm the correct container engine is enabled:

```
-[nf-core/rnaseq] Singularity enabled
```

If Nextflow detects the wrong container engine or none at all, verify your container scope settings.

### Job submission

Monitor job submissions to your cluster scheduler:

```bash
# For SLURM
squeue -u $USER

# For SGE
qstat -u $USER

# For PBS
qstat -u $USER
```

Verify that:

- Jobs submit to the correct queue/partition
- Resource requests match your configuration
- Jobs include any custom cluster options you specified

### Job execution

Watch for jobs to complete successfully:

```bash
tail -f .nextflow.log
```

Check the Nextflow log for errors related to:

- Resource allocation
- Container image pulling
- File system access
- Module loading (if using `beforeScript`)

### Work directory

Verify that intermediate files appear in the expected location:

```bash
ls work/
```

If using scratch directories, confirm files write to the configured scratch path rather than your home directory.

## Test different resource requirements

Test with pipelines that stress different aspects of your configuration:

### Small resource pipeline

Test basic functionality with minimal resources:

```bash
nextflow run nf-core/fetchngs \
  -profile big_university,test \
  --custom_config_base /path/to/your/configs
```

### Large resource pipeline

Test resource limit handling:

```bash
nextflow run nf-core/sarek \
  -profile big_university,test \
  --custom_config_base /path/to/your/configs
```

### Multiple sub-profiles

If you defined sub-profiles, test each one:

```bash
nextflow run nf-core/rnaseq \
  -profile big_university,bigmem,test \
  --custom_config_base /path/to/your/configs
```

## Test retry behavior

Verify that job retries work correctly by checking the Nextflow log for retry messages. The `maxRetries` setting should allow jobs to resubmit with increased resources when they fail.

## Troubleshooting

### Profile not loading

**Problem:** Your profile does not appear in the pipeline startup output

**Solution:** Verify your profile configuration and loading:

- Check for syntax errors in `conf/<cluster_name>.config`
- Verify your profile is listed in `nfcore_custom.config`
- Ensure you specified the correct profile name with `-profile`

### Wrong container engine

**Problem:** Nextflow detects the wrong container engine or none at all

**Solution:** Check your container configuration:

- Only enable one container engine in your configuration
- Verify the container executable is available on your system with `singularity --version` or `docker --version`
- Confirm the container scope name matches the engine name (for example, `singularity`, not `apptainer`)

### Jobs not submitting

**Problem:** Jobs do not submit to your scheduler

**Solution:** Verify your scheduler configuration:

- Check that the `executor` value matches your scheduler (for example, `slurm`, `sge`, `pbs`)
- Ensure Nextflow can access the scheduler commands (`sbatch`, `qsub`, etc.) by running them directly
- Review `clusterOptions` for syntax errors
- Verify you have permissions to submit jobs to the configured queue

### Resource limit errors

**Problem:** Jobs fail with out-of-memory or resource allocation errors

**Solution:** Adjust your resource configuration:

- Verify `resourceLimits` values match your cluster's maximum node specifications
- Check that queue assignments can accommodate the requested resources
- Test with the `test` profile to use minimal resources and isolate the issue
- Review queue-specific limits and ensure your configuration respects them

## Next steps

After successful testing, you are ready to submit your profile to the nf-core/configs repository.

<!-- TODO: Link to page about making PRs -->

