---
title: Configuration file components
subtitle: Structuring your institutional profile configuration
shortTitle: Configuration
weight: 2
---

The main configuration file (`conf/<cluster_name>.config`) defines all cluster-specific settings using Nextflow's configuration syntax.
This file contains multiple configuration scopes that control different aspects of pipeline execution.

## params scope

The `params` scope contains nf-core-specific metadata and shared resource paths.
These parameters appear in pipeline execution logs and help users identify which profile loaded.

```groovy
params {
  config_profile_description = 'Big University HPC cluster profile provided by nf-core/configs'
  config_profile_contact = 'Jane Smith (@jsmith)'
  config_profile_url = 'https://hpc.biguniversity.edu'
}
```

### Required parameters

Include these parameters in every institutional profile:

- `config_profile_description`: Brief description of your profile and institution
- `config_profile_contact`: Name and GitHub handle of the profile maintainer
- `config_profile_url`: URL to your institution's HPC documentation or homepage

### Optional parameters

Include these parameters when applicable:

- `igenomes_base`: Path to a local iGenomes reference genome directory
- `max_memory`: Maximum memory available on the largest node (for example, `2.TB`, `768.GB`)
- `max_cpus`: Maximum CPU cores available on the largest node (for example, `128`, `64`)
- `max_time`: Maximum walltime for jobs (for example, `720.h`, `168.h`)

:::note
The `max_memory`, `max_cpus`, and `max_time` parameters are legacy parameters that some older pipelines still reference.
Newer pipelines use `resourceLimits` in the `process` scope instead.
Include both for maximum compatibility.
:::

### Example params scope

```groovy
params {
  config_profile_description = 'Big University HPC cluster profile provided by nf-core/configs'
  config_profile_contact = 'Jane Smith (@jsmith)'
  config_profile_url = 'https://hpc.biguniversity.edu'
  igenomes_base = '/shared/data/igenomes/'

  // Legacy parameters for compatibility
  max_memory = 2.TB
  max_cpus = 128
  max_time = 720.h
}
```

## process scope

The `process` scope defines resource limits, executor settings, and job submission details.
This is where you configure how Nextflow submits jobs to your cluster scheduler.

### Resource limits

Define maximum resources using the `resourceLimits` directive:

```groovy
process {
  resourceLimits = [
    cpus: 128,
    memory: 2.TB,
    time: 720.h
  ]
}
```

These limits prevent pipelines from requesting more resources than your cluster can provide.
When a process requests resources beyond these limits, Nextflow automatically caps the request at the maximum value.

### Executor and queue

Specify the job scheduler and default queue:

```groovy
process {
  executor = 'slurm'
  queue = 'general'
}
```

Common executor values:

- `slurm`: SLURM scheduler
- `sge`: Sun Grid Engine
- `pbs`: PBS/Torque
- `lsf`: IBM LSF
- `local`: Local execution (no scheduler)

### Retry strategy

Configure automatic job retry with exponentially increasing resources:

```groovy
process {
  maxRetries = 2
}
```

When a job fails due to resource limits, Nextflow automatically resubmits it with more resources.
Most profiles set `maxRetries` to 2, allowing jobs to attempt execution with progressively higher resource allocations.

### Dynamic queue selection

Use Groovy expressions to route jobs to appropriate queues based on requirements:

```groovy
process {
  queue = { task.time <= 1.h ? 'short' : task.time <= 24.h ? 'medium' : 'long' }
}
```

This example routes jobs to different queues based on requested runtime:

- Jobs requesting ≤ 1 hour go to the `short` queue
- Jobs requesting ≤ 24 hours go to the `medium` queue
- Jobs requesting > 24 hours go to the `long` queue

You can also route based on CPU or memory requirements:

```groovy
process {
  queue = { task.cpus > 64 ? 'highmem' : 'general' }
}
```

### Additional process directives

Other useful directives for the process scope:

- `clusterOptions`: Add scheduler-specific options (for example, `clusterOptions = '--account=projectA'`)
- `beforeScript`: Load modules or set environment variables before job execution
- `scratch`: Specify temporary directory for intermediate files

### Example process scope

```groovy
process {
  resourceLimits = [
    cpus: 128,
    memory: 2.TB,
    time: 720.h
  ]
  executor = 'slurm'
  queue = { task.time <= 24.h ? 'general' : 'long' }
  maxRetries = 2
  clusterOptions = '--account=nf_core_users'
}
```

## executor scope

The `executor` scope controls how Nextflow submits jobs to your cluster scheduler.
These settings prevent overwhelming small clusters with too many simultaneous submissions.

```groovy
executor {
  queueSize = 8
  submitRateLimit = '10 sec'
}
```

### Executor parameters

- `queueSize`: Maximum number of jobs Nextflow can submit simultaneously
- `submitRateLimit`: Minimum time between job submissions

Adjust these values based on your cluster size and scheduler policies.
Smaller clusters benefit from lower values to prevent scheduler overload.

### Example executor scope

```groovy
executor {
  queueSize = 10
  submitRateLimit = '5 sec'
}
```

## Container scopes

Container scopes configure container engine settings.
Choose the scope that matches your cluster's available container engine.

### Singularity configuration

Most HPC clusters use Singularity (now called Apptainer):

```groovy
singularity {
  enabled = true
  autoMounts = true
  cacheDir = '/shared/containers/singularity'
}
```

Singularity parameters:

- `enabled`: Enables Singularity for all processes
- `autoMounts`: Automatically mounts common system paths
- `cacheDir`: Shared directory for storing container images

### Docker configuration

For clusters with Docker support:

```groovy
docker {
  enabled = true
  runOptions = '-u $(id -u):$(id -g)'
}
```

Docker parameters:

- `enabled`: Enables Docker for all processes
- `runOptions`: Additional options for `docker run` commands

### Podman configuration

For clusters using Podman:

```groovy
podman {
  enabled = true
}
```

### Conda configuration

For clusters that prefer Conda environments:

```groovy
conda {
  enabled = true
  cacheDir = '/shared/conda/envs'
}
```

:::tip
Only enable one container engine per profile.
Users can override the container engine with their own profile if needed.
:::

## profiles scope

The `profiles` scope defines internal sub-profiles for cluster variants.
Use this when your institution has multiple clusters or partitions with different resource limits.

```groovy
profiles {
  bigmem {
    process {
      resourceLimits = [
        cpus: 64,
        memory: 4.TB,
        time: 168.h
      ]
      queue = 'bigmem'
    }
  }

  gpu {
    process {
      resourceLimits = [
        cpus: 32,
        memory: 512.GB,
        time: 72.h
      ]
      queue = 'gpu'
      clusterOptions = '--gres=gpu:1'
    }
  }
}
```

Users can access sub-profiles by combining them with your main profile:

```bash
nextflow run nf-core/rnaseq -profile big_university,bigmem
```

:::warning
Internal sub-profiles do not inherit settings from the main configuration.
You must redefine all necessary configuration scopes within each sub-profile.
:::

## Complete configuration example

Here is a complete example combining all scopes:

```groovy
// Params scope - metadata and paths
params {
  config_profile_description = 'Big University HPC cluster profile provided by nf-core/configs'
  config_profile_contact = 'Jane Smith (@jsmith)'
  config_profile_url = 'https://hpc.biguniversity.edu'
  igenomes_base = '/shared/data/igenomes/'

  // Legacy parameters
  max_memory = 2.TB
  max_cpus = 128
  max_time = 720.h
}

// Process scope - resource limits and scheduler settings
process {
  resourceLimits = [
    cpus: 128,
    memory: 2.TB,
    time: 720.h
  ]
  executor = 'slurm'
  queue = { task.time <= 24.h ? 'general' : 'long' }
  maxRetries = 2
  beforeScript = 'module load nextflow'
}

// Executor scope - submission controls
executor {
  queueSize = 10
  submitRateLimit = '5 sec'
}

// Singularity scope - container settings
singularity {
  enabled = true
  autoMounts = true
  cacheDir = '/shared/containers/singularity'
}

// Internal profiles for cluster variants
profiles {
  bigmem {
    process {
      resourceLimits = [
        cpus: 64,
        memory: 4.TB,
        time: 168.h
      ]
      queue = 'bigmem'
    }
  }
}
```

## Next steps

After creating your configuration file, document your profile to help users at your institution.

See [Documentation requirements](./documentation.md) to learn what to include in your documentation file.
