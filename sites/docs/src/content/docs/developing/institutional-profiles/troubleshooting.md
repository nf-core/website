---
title: Troubleshooting
subtitle: Troubleshoot institutional profiles
shortTitle: Troubleshooting
---

This page covers common issues you may encounter when developing, running, or testing institutional profiles.

## Testing

### Profile not loading

**Problem:** Your profile does not appear in the pipeline startup output

**Solution:** Verify your profile configuration and loading:

- Check for syntax errors in `conf/<cluster-name>.config`
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
