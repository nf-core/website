---
title: Documentation requirements
subtitle: Learn what to include in your institutional profile documentation
shortTitle: Documentation
weight: 5
---

Each institutional profile requires a documentation file (`docs/<cluster_name>.md`) that explains how to use the profile and provides cluster-specific information. Clear documentation helps users at your institution get started quickly and troubleshoot common issues.

## Required sections

Your documentation file should include these sections at minimum:

### Profile description

Start with a brief introduction that identifies the cluster and institution:

```markdown
# Big University HPC Configuration

Configuration profile for the Big University HPC cluster.

This profile is maintained by the Research Computing team at Big University.
Contact us at [hpc-support@biguniversity.edu](mailto:hpc-support@biguniversity.edu)
for questions or issues.
```

### Usage instructions

Explain how users at your institution should use the profile:

````markdown
## Using the profile

To use this profile, specify it when running any nf-core pipeline:

\```bash
nextflow run nf-core/<pipeline_name> -profile big_university
\```

This will automatically apply all cluster-specific settings including executor,
resource limits, and container configurations.
````

Include any sub-profiles if you defined them:

````markdown
### Sub-profiles

The following sub-profiles are available for specialized resources:

- `big_university,bigmem`: Use the high-memory partition (4TB RAM nodes)
- `big_university,gpu`: Use the GPU partition with CUDA support

Example:

\```bash
nextflow run nf-core/rnaseq -profile big_university,bigmem
\```
````

### Prerequisites and setup

Document any setup steps users must complete before running pipelines:

```markdown
## Prerequisites

Before running pipelines with this profile:

1. Ensure you have access to the Big University HPC cluster
1. Load the Nextflow module: `module load nextflow`
1. Verify Singularity is available: `singularity --version`

### First-time setup

The first time you run a pipeline, Nextflow will download container images to
`/shared/containers/singularity`. This may take several minutes depending on
the pipeline. Subsequent runs will reuse these cached images.
```

### Cluster information

Provide relevant details about your cluster:

```markdown
## Cluster details

- **Scheduler**: SLURM
- **Default queue**: `general`
- **Container engine**: Singularity
- **Resource limits**:
  - Maximum CPUs: 128 cores
  - Maximum memory: 2TB
  - Maximum walltime: 720 hours
```

## Optional sections

Include these sections when applicable:

### Queue information

If your cluster has multiple queues with different characteristics:

```markdown
## Available queues

| Queue   | Max CPUs | Max Memory | Max Time | Purpose          |
| ------- | -------- | ---------- | -------- | ---------------- |
| general | 128      | 2TB        | 720h     | Standard jobs    |
| short   | 64       | 512GB      | 1h       | Quick jobs       |
| long    | 64       | 1TB        | 2160h    | Extended jobs    |
| bigmem  | 64       | 4TB        | 168h     | Memory-intensive |

The profile automatically routes jobs to appropriate queues based on requirements.
```

### Reference genomes

If you provide local reference genome paths:

```markdown
## Reference genomes

This profile includes paths to local iGenomes reference genomes at
`/shared/data/igenomes/`. Pipelines will automatically use local copies
instead of downloading references, saving time and bandwidth.
```

### Known issues and limitations

Document any cluster-specific issues or limitations:

```markdown
## Known issues

- Jobs requesting more than 24 hours may experience longer queue times
- The `/home` filesystem has a 100GB quota. Use `/scratch` for large datasets
- Container builds are not supported. Contact support if you need custom containers
```

### Troubleshooting

Include common problems and solutions:

````markdown
## Troubleshooting

### Out of memory errors

If jobs fail with out-of-memory errors, try the `bigmem` sub-profile:

\```bash
nextflow run nf-core/<pipeline_name> -profile big_university,bigmem
\```

### Quota exceeded errors

Move your work directory to `/scratch` if you exceed `/home` quota:

\```bash
cd /scratch/$USER
nextflow run nf-core/<pipeline_name> -profile big_university
\```
````

### Contact information

Provide support contacts for users:

```markdown
## Support

For issues related to this profile or the Big University HPC cluster:

- Email: [hpc-support@biguniversity.edu](mailto:hpc-support@biguniversity.edu)
- Documentation: https://hpc.biguniversity.edu/docs
- Submit tickets: https://hpc.biguniversity.edu/support

For issues with nf-core pipelines, see the [nf-core website](https://nf-co.re).
```

## Complete documentation example

Here is a complete example combining all sections:

````markdown
# Big University HPC Configuration

Configuration profile for the Big University HPC cluster.

This profile is maintained by the Research Computing team at Big University.
Contact us at [hpc-support@biguniversity.edu](mailto:hpc-support@biguniversity.edu)
for questions or issues.

## Using the profile

To use this profile, specify it when running any nf-core pipeline:

\```bash
nextflow run nf-core/<pipeline_name> -profile big_university
\```

This will automatically apply all cluster-specific settings including executor,
resource limits, and container configurations.

## Prerequisites

Before running pipelines with this profile:

1. Ensure you have access to the Big University HPC cluster
1. Load the Nextflow module: `module load nextflow`
1. Verify Singularity is available: `singularity --version`

## Cluster details

- **Scheduler**: SLURM
- **Default queue**: `general`
- **Container engine**: Singularity
- **Resource limits**:
  - Maximum CPUs: 128 cores
  - Maximum memory: 2TB
  - Maximum walltime: 720 hours

## Reference genomes

This profile includes paths to local iGenomes reference genomes at
`/shared/data/igenomes/`. Pipelines will automatically use local copies
instead of downloading references.

## Troubleshooting

### Out of memory errors

If jobs fail with out-of-memory errors, contact support about accessing
the high-memory partition.

### Quota exceeded errors

Move your work directory to `/scratch` if you exceed `/home` quota:

\```bash
cd /scratch/$USER
nextflow run nf-core/<pipeline_name> -profile big_university
\```

## Support

For issues related to this profile or the Big University HPC cluster:

- Email: [hpc-support@biguniversity.edu](mailto:hpc-support@biguniversity.edu)
- Documentation: https://hpc.biguniversity.edu/docs

For issues with nf-core pipelines, see the [nf-core website](https://nf-co.re).
````

## Next steps

After creating your documentation file, test your profile to verify all settings work correctly.

See [Testing profiles](./testing.md) to learn how to test your institutional profile.
