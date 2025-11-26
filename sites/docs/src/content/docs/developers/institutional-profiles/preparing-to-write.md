---
title: Preparing to write a profile
subtitle: Gather the information you need before creating an institutional profile
shortTitle: Preparing to write
weight: 2
---

Before creating an institutional profile, gather critical information about your cluster and confirm that you have permission to share configuration details publicly.

## Check for existing profiles

Before starting, verify that a profile does not already exist for your institution:

1. Browse the [nf-core/configs repository](https://github.com/nf-core/configs/tree/master/conf) to see existing profiles
1. Check the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation) for a complete list
1. Search for your institution name or cluster name in the repository

If a profile already exists, consider contributing improvements rather than creating a duplicate.

## Determine the appropriate scope

Choose the right level for your institutional profile:

- **Single cluster**: Create one profile for a specific cluster (for example, `bigcluster` or `hpc_system`)
- **Multiple clusters**: Create sub-profiles within a parent institution profile when your organization has multiple distinct clusters with different configurations
- **Umbrella institution**: Consider whether a single profile can serve multiple related systems or if separate profiles provide better clarity

The profile name should clearly identify the institution or cluster. Avoid generic names like `university` or `hpc` that could apply to many institutions.

## Confirm public sharing permissions

Institutional profiles are publicly visible in the nf-core/configs repository. Before proceeding:

1. Verify with your system administrators that cluster configuration details can be shared publicly
1. Confirm that resource limits, queue names, and file paths do not expose sensitive information
1. Ensure that you have permission to represent your institution in the nf-core community

Most cluster configurations contain non-sensitive information, but institutional policies vary.

## Gather cluster information

Collect the following information about your cluster. Contact your system administrators if you do not have access to these details.

### Scheduler and executor

Identify the job scheduling system your cluster uses:

- SLURM
- SGE (Sun Grid Engine)
- PBS/Torque
- LSF
- Moab
- HTCondor
- Local execution (for small systems)

This determines the `executor` value in your configuration and affects how Nextflow submits jobs.

### Container support

Determine which container engines are available on your cluster:

- **Singularity/Apptainer**: Most common on HPC systems
- **Docker**: Common on cloud platforms and local systems
- **Podman**: Alternative to Docker with rootless support
- **Conda**: Software environment management (not containerization)

Note whether users need special permissions or module loads to access these tools.

### Resource maximums

Document the largest node specifications in your cluster:

- Maximum memory per node (for example, `2.TB`, `768.GB`)
- Maximum CPU cores per node (for example, `128`, `64`)
- Maximum walltime per job (for example, `720.h`, `168.h`)

These values define the upper limits in your profile's `resourceLimits` configuration. Base these on the largest single node specifications, not aggregate cluster totals.

### Queues and partitions

List all available queues or partitions and their specifications:

- Queue/partition names
- Resource limits for each queue (CPUs, memory, walltime)
- Access restrictions (for example, specific user groups)
- Purpose or typical use cases

This information helps you configure default queues and set up dynamic queue selection based on job requirements.

### Shared resources

Note paths to shared resources that benefit all users:

- Reference genome directories (iGenomes paths if available)
- Container/Singularity image cache directories
- Software module locations
- Shared scratch or temporary directories

Sharing these paths reduces redundant downloads and storage usage across users.

### Scratch directory requirements

Determine if your cluster requires per-node temporary storage:

- Does the scheduler provide a per-job temporary directory?
- Do processes need to write to local node storage rather than shared filesystems?
- What is the typical scratch directory path (for example, `/tmp`, `/scratch`, `$TMPDIR`)?

Some clusters have strict quotas on shared filesystems, making scratch directories essential for intermediate files.

## Additional considerations

Consider these additional factors that may affect your profile:

- **Module system**: Does your cluster use Environment Modules or Lmod? Note which modules users must load before running Nextflow
- **Network configuration**: Are there specific network zones or storage systems that affect job execution?
- **Billing or accounting**: Does your cluster require specific account codes or project identifiers in job submissions?
- **Submit rate limits**: Does your scheduler have restrictions on how quickly you can submit jobs?

## Next steps

Once you have gathered this information, you are ready to create the required files for your institutional profile.

See [File structure requirements](./file-structure.md) to learn which files to create or modify.
