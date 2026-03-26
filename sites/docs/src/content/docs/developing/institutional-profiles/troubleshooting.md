---
title: Troubleshooting
subtitle: Troubleshoot institutional profiles
shortTitle: Troubleshooting
---

This page covers common issues you may encounter when developing, running, or testing institutional profiles.

## Testing

### Profile not loading

A profile may not appear in the pipeline startup output.

This issue occurs when the profile configuration is invalid or not correctly registered. To resolve this issue:

- Check for syntax errors in `conf/<cluster-name>.config`
- Verify your profile is listed in `nfcore_custom.config`
- Ensure you specified the correct profile name with `-profile`

### Wrong container engine detected

Nextflow may detect the wrong container engine or none at all.

This issue occurs when multiple container engines are enabled or the container scope name does not match the engine. To resolve this issue:

- Only enable one container engine in your configuration
- Verify the container executable is available on your system with `singularity --version` or `docker --version`
- Confirm the container scope name matches the engine name (for example, `singularity`, not `apptainer`)

### Jobs not submitting to scheduler

Jobs may not submit to your scheduler.

This issue occurs when the executor is misconfigured or Nextflow cannot access the scheduler commands. To resolve this issue:

- Check that the `executor` value matches your scheduler (for example, `slurm`, `sge`, `pbs`)
- Ensure Nextflow can access the scheduler commands (`sbatch`, `qsub`, etc.) by running them directly
- Review `clusterOptions` for syntax errors
- Verify you have permissions to submit jobs to the configured queue

### Resource limit errors

Jobs may fail with out-of-memory or resource allocation errors.

This issue occurs when `resourceLimits` values exceed what the cluster can accommodate, or when queue assignments do not match the requested resources. To resolve this issue:

- Verify `resourceLimits` values match your cluster's maximum node specifications
- Check that queue assignments can accommodate the requested resources
- Test with the `test` profile to use minimal resources and isolate the issue
- Review queue-specific limits and ensure your configuration respects them
