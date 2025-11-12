---
title: "Managing Nextflow work directory growth"
subtitle: "A guide for efficient storage utilization"
---

The management of intermediate files generated during Nextflow pipeline execution is a challenge for some workflows. As pipelines increase in complexity and scale, work directories can rapidly consume available storage, potentially leading to pipeline failures. This guide summarizes strategies for managing work directory growth while maintaining pipeline reproducibility and debugging capabilities.

## Work directory accumulation

The Nextflow work directory serves as an important component of the execution model, providing caching capabilities and resume functionality. During pipeline execution, Nextflow creates unique subdirectories for each task (e.g., `work/3f/70944c7a549b6221e1ccc7b4b21b62`) containing:

- Symbolic links to input files
- Intermediate output files
- Command scripts and execution logs
- Temporary files created during tasks

In production environments processing large-scale genomic datasets, work directories can expand rapidly.

## Storage management options

The following sections describe ways to manage the storage of files.

### Selective post-execution cleanup

Nextflow's built-in `clean` command enables targeted removal of work directories. The following command preserves work directories from the current execution while removing directories from previous executions:

```bash
nextflow clean -f -before $(nextflow log -q | tail -n 1)
```

Command components:

- `nextflow log -q{:bash}`: Returns a list of run names
- `tail -n 1{:bash}`: Isolates the most recent execution identifier
- `-before{:bash}`: Specifies cleanup of executions preceding the specified run
- `-f{:bash}`: Executes deletion without confirmation

For verification, perform a dry run using the `-n{:bash}` option:

```bash
nextflow clean -n -before $(nextflow log -q | tail -n 1)
```

### Automated cleanup configuration

Nextflow supports automatic work directory cleanup upon successful pipeline completion through configuration directives:

```groovy title="nextflow.config"
cleanup = true
```

:::note
Enabling automatic cleanup prevents the use of resume functionality for the affected pipeline execution. This configuration suits production pipelines where output reproducibility is assured and resume capability isn't required.
:::

### Scratch directory implementation

The scratch directive enables process execution in temporary directories, typically local to compute nodes, with selective output staging to the work directory:

```groovy
process SEQUENCE_ALIGNMENT {
    scratch true

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.bam")

    script:
    """
    bwa mem reference.fa ${reads} > aligned.sam
    samtools sort aligned.sam > sorted.bam
    samtools index sorted.bam
    # Only sorted.bam transfers to work directory
    """
}
```

:::tip
This configuration is particularly beneficial in HPC environments where it reduces network filesystem overhead.
:::

### Dynamic intermediate file management with nf-boost

The [`nf-boost`](https://registry.nextflow.io/plugins/nf-boost) plugin implements intelligent cleanup mechanisms that remove intermediate files during pipeline execution as they become unnecessary:

```groovy title="nextflow.config"
plugins {
    id 'nf-boost'
}

boost {
    cleanup = true
    cleanupInterval = '180s'  // Cleanup evaluation interval
}
```

See [nf-boost](https://github.com/bentsherman/nf-boost) for more information.

### Pipeline optimization strategies

Minimize intermediate file generation through process optimization:

```groovy
process OPTIMIZED_ANALYSIS {
    input:
    path input_data

    output:
    path "final_results.txt"

    script:
    """
    # Utilize pipe operations to avoid intermediate files
    initial_process ${input_data} | \\
    intermediate_transform | \\
    final_analysis > final_results.txt

    # Implement named pipes for file-dependent tools
    mkfifo temp_pipe
    producer_process ${input_data} > temp_pipe &
    consumer_process temp_pipe > final_results.txt
    rm temp_pipe

    # Use process substitution to avoid intermediate files
    paste <( cat ${input_data} ) <( cat ${input_data} ) > final_results.txt
    """
}
```

## Recommendations

Effective management of Nextflow work directories requires a tailored approach. The Nextflow `clean` command provides essential functionality for storage recovery and work directory maintenance, though implementation must balance storage optimization with requirements for pipeline resumption and debugging.

For development environments, combine reduced test datasets with manual cleanup for optimal flexibility. Production deployments benefit from automatic cleanup or dynamic solutions like `nf-boost`. HPC installations should leverage scratch directories to minimize shared storage impact.

You should establish storage management policies incorporating:

- Regular maintenance schedules
- Environment-specific configuration profiles
- Capacity planning procedures
- Documentation of cleanup strategies

Through systematic implementation of these strategies, you can maintain efficient pipeline operations while preventing storage exhaustion, ensuring sustainable computational workflow execution at scale.
