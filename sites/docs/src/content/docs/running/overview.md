---
title: Overview
subtitle: Learn how to run and configure nf-core pipelines
shortTitle: Overview
weight: 1
---

<!-- TODO: Verify all pages in the running folder are mentioned, sanity check content accuracy, and validate all links -->

nf-core pipelines are designed to be accessible, flexible, and reproducible across different computing environments.
Whether you're running a pipeline on a laptop, HPC cluster, or cloud infrastructure, nf-core provides consistent command structures and comprehensive configuration options to meet your needs.

This section covers everything you need to run nf-core pipelines effectively, from basic execution commands to advanced configuration and specialised scenarios.

## Getting started

All nf-core pipelines follow a consistent command structure and execution pattern.
Start here to learn the essential commands for running any nf-core pipeline with your own data.

- **[Running pipelines](/running/run-a-pipeline):** Essential commands and patterns for running nf-core pipelines, including testing with sample data, using parameter files, and resuming failed runs

## Configuration

nf-core pipelines can be configured to work with different execution environments, resource requirements, and infrastructure constraints.

- **[Configuration options](/running/configuration/configuration-options):** Detailed guidance on configuring options
- **[System requirements](/running/configuration/nextflow-for-your-system):** Guidance on configuring pipelines to match your system's capabilities, including resource allocation, executors, and tool arguments

## Reference data

Many nf-core pipelines require reference genomes and annotation files.
Learn how to access and manage reference data efficiently.

- **[Reference genomes](/running/reference-genomes):** Approaches for managing reference genomes, including AWS iGenomes, custom genome files, and Refgenie

## Running pipelines offline

For systems without internet access, nf-core provides solutions for preparing and transferring all required components.

- **[Running pipelines offline](/running/run-pipelines-offline):** Guidance on preparing and running nf-core pipelines on systems without internet access, including transferring pipeline code, containers, and reference data

## Advanced topics

For specialised computing environments or resource management requirements, these guides address specific challenges in pipeline execution.

- **[Google Colab](/running/advanced-topics/google-colab):** Guidance on running nf-core pipelines using Google Colab's cloud resources, addressing limitations in local computing environments
- **[Managing work directory growth](/running/advanced-topics/managing_work_directory_growth):** Strategies for managing intermediate files and work directory storage during pipeline execution
