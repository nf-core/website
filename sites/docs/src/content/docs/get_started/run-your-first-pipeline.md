---
title: Run your first pipeline
subtitle: Run the nf-core/demo pipeline
shortTitle: Your first pipeline
weight: 3
---

With Nextflow and your software dependency manager installed, you are ready to execute your first nf-core pipeline.

## nf-core/demo

The `nf-core/demo` pipeline is a lightweight demonstration pipeline designed to introduce the nf-core framework. It performs basic quality control and processing on sequencing data, showcasing the structure and best practices followed by all nf-core pipelines. This pipeline is ideal for verifying your environment is correctly configured before working with production pipelines.

## Running the pipeline

Execute the pipeline using a test dataset:

```bash
nextflow run nf-core/demo -r 1.0 -profile test,<software_config>
```

Replace `<software_config>` with your software dependency manager: `docker`, `singularity`, `podman`, or `conda`.

For example, to run with Docker:

```bash
nextflow run nf-core/demo -r 1.0 -profile test,docker
```

The `test` profile uses a small dataset bundled with the pipeline, allowing you to verify that the pipeline executes correctly without requiring your own data.

## Understanding the command

Each component of the command serves a specific purpose:

- `nextflow run`: The Nextflow command to execute a pipeline
- `nf-core/demo`: The pipeline identifier, composed of the GitHub organization (`nf-core`) and repository name (`demo`)
- `-r 1.0`: Specifies the pipeline version to use
- `-profile test,<software_config>`: Configuration profiles that define pipeline behavior:
  - `test`: Configures the pipeline to use bundled test data and parameters
  - `<software_config>`: Specifies the software dependency management system

When you execute the command, Nextflow will automatically download the pipeline from GitHub and manage its execution.

## Viewing results

Pipeline outputs will be written to the `results` directory in your current working location. The results include quality control reports, processed data, and execution logs.
