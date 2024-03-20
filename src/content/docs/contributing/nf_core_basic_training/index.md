---
title: Basic training to create an nf-core pipeline
subtitle: A guide to create Nextflow pipelines using nf-core tools
---

## Scope

- How do I create a pipeline using nf-core tools?
- How do I incorporate modules from nf-core modules?
- How can I use custom code in my pipeline?

:::note

### Learning objectives

- The learner will create a simple pipeline using the nf-core template.
- The learner will identify key files in the pipeline.
- The learner will lint their pipeline code to identify work to be done.
- The learner will incorporate modules from nf-core/modules into their pipeline.
- The learner will add custom code as a local module into their pipeline.
- The learner will build an nf-core schema to describe and validate pipeline parameters.

:::

This training course aims to demonstrate how to build an nf-core pipeline using the nf-core pipeline template and nf-core modules and subworkflows as well as custom, local modules. Be aware that we are not going to explain any fundamental Nextflow concepts, as such we advise anyone taking this course to have completed the [Basic Nextflow Training Workshop](https://training.nextflow.io/).

```md
During this course we are going to build a Simple RNA-Seq workflow.
This workflow is by no means ment to be a useful bioinformatics workflow,
but should only teach the objectives of the course, so please,
**DO NOT use this workflow to analyse RNA sequencing data**!
```

### Follow the training videos

This training can be followed either based on this documentation alone, or via a training video hosted on youtube. You can find the youtube video in the Youtube playlist below:

(no such video yet)

## Overview

### Layout of the pipeline

The course is going to build an (totally unscientific and useless) RNA seq pipeline that does the following:

1. Indexing of a transcriptome file
2. Quality control
3. Quantification of transcripts
4. [whatever the custom script does]
5. Generation of a MultiQC report

### Let's get started

1. **Setting up the gitpod environment for the course**

The course is using gitpod in order to avoid the time expense for downloading and installing tools and data. [Learn how to setup the GitPod environment](/docs/contributing/nf_core_basic_training/gitpod_environment.md)

2. **Creating a new nf-core pipeline from the nf-core template**

   a) generate the pipeline with `nf-core create`

   b) The template git repository

   c) Running the pipeline using the test profile

   d) Linting the pipeline

   [These steps are described in the section "Generate a pipeline with nf-core tools"](/docs/contributing/nf_core_basic_training/nf_core_create_tool.md)

3. **Exploring the nf-core template files**

   The template contains a range of important files and directories. Check them out in the [walk-through of the template files](/docs/contributing/nf_core_basic_training/template_walk_through.md)

4. **Building a nf-core pipeline using the template**

   a) [Adding a nf-core module to your pipeline](/docs/contributing/nf_core_basic_training/add_nf_core_module.md)

   b) [Adding a local custom module to your pipeline](/docs/contributing/nf_core_basic_training/add_custom_module.md)

   c) [Working with Nextflow schema](/docs/contributing/nf_core_basic_training/nf_schema.md)

   d) [Linting your modules](/docs/contributing/nf_core_basic_training/linting_modules.md)

:::note

### Key points

- `nf-core create <pipeline>` creates a pipeline from the nf-core template.
- `nf-core lint` lints the pipeline code for things that must be completed.
- `nf-core modules list local` lists modules currently installed into your pipeline.
- `nf-core modules list remote` lists modules available to install into your pipeline.
- `nf-core modules install <tool/subtool>` installs the tool module into your pipeline.
- `nf-core modules create` creates a module locally to add custom code into your pipeline.
- `nf-core modules lint --all` lints your module code for things that must be completed.
- `nf-core schema build` opens an interface to allow you to describe your pipeline parameters and set default values, and which values are valid.

:::
