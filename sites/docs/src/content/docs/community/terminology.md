---
title: Terminology
subtitle: Specification of terms used in the nf-core community
shortTitle: Terminology
---

Nextflow DSL2 lets you write pipelines at different levels of granularity. This page defines the terminology for DSL2 components in nf-core.

## Domain-Specific Language (DSL)

A domain-specific language (DSL) is a programming language developed for a specific application. Nextflow uses DSL2, the latest version, which supports scalable and modular data analysis pipelines.

## Module

A `process` that you can use across different pipelines. Modules are atomic and cannot be split further. For example, a module file contains the process definition for a single tool such as `FastQC`. nf-core module files are available in the [`modules/`](https://github.com/nf-core/modules/tree/master/modules) directory of nf-core/modules with documentation and tests.

## Subworkflow

A chain of multiple modules that provides higher-level functionality in a pipeline. For example, a subworkflow runs multiple QC tools with FastQ files as input. You should ship subworkflows with the pipeline implementation and share them between pipelines when needed. nf-core subworkflow files are available in the [`subworkflow/`](https://github.com/nf-core/modules/tree/master/subworkflows) directory of nf-core/modules with documentation and tests.

## Component

Component is our umbrella term for both modules and subworkflows, because they share a lot of the same characteristics.

## Workflow

An end-to-end pipeline where one or more inputs produce a series of outputs. You can implement workflows using a large monolithic script or a combination of DSL2 modules and subworkflows. nf-core pipelines can have multiple workflows. For example, [nf-core/viralrecon](https://github.com/nf-core/viralrecon/tree/master/workflows) processes different data types for the same purpose.
