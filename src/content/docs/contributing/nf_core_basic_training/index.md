---
title: Basic training to create an nf-core pipeline
subtitle: A guide to create Nextflow pipelines using nf-core tools
---

# Scope

- How do I create a pipeline using nf-core tools?
- How do I incorporate modules from nf-core modules?
- How can I use custom code in my pipeline?

:::note

## Learning objectives

- The learner will create a simple pipeline using the nf-core template.
- The learner will identify key files in the pipeline.
- The learner will lint their pipeline code to identify work to be done.
- The learner will incorporate modules from nf-core/modules into their pipeline.
- The learner will add custom code as a local module into their pipeline.
- The learner will build an nf-core schema to describe and validate pipeline parameters.

:::

This training course aims to demonstrate how to build an nf-core pipeline using the nf-core pipeline template and nf-core modules and subworkflows as well as custom, local modules. Be aware that we are not going to explain any fundamental Nextflow concepts, as such we advise anyone taking this course to have completed the [Basic Nextflow Training Workshop](https://training.nextflow.io/).

# Let's get started

## Follow the training videos

This training can be followed either based on this documentation alone, or via a training video hosted on youtube. You can find the youtube video in the Youtube playlist below:

(no such video yet)

## Training sessions

<p class="text-left">
  <a href="gitpod_environment/" class="btn btn-lg btn-success" target="_blank">
    Chapter 1 - setup the GitPod environment
  </a>
</p>

The course is using gitpod in order to avoid the time expense for downloading and installing tools and data.

<p class="text-left">
  <a href="nf_core_create_tool/" class="btn btn-lg btn-success" target="_blank">
    Chapter 2 - Creating a new nf-core pipeline from the nf-core template
  </a>
</p>

This includes the following sub-sections:

a) generate the pipeline with `nf-core create`

b) The template git repository

c) Running the pipeline using the test profile

d) Linting the pipeline

<p class="text-left">
  <a href="nf_core_create_tool/" class="btn btn-lg btn-success" target="_blank">
    Chapter 3 - Exploring the nf-core template files
  </a>
</p>

The template contains a range of important files and directories. This chapter is giving explanations to all the template content important to beginners.

<p class="text-left">
  <a href="add_nf_core_module/" class="btn btn-lg btn-success" target="_blank">
    Chapter 4 - Adding a nf-core module to your pipeline
  </a>
</p>
<p class="text-left">
  <a href="add_custom_module/" class="btn btn-lg btn-success" target="_blank">
    Chapter 5 - Adding a local custom module to your pipeline
  </a>
</p>
<p class="text-left">
  <a href="nf_schema/" class="btn btn-lg btn-success" target="_blank">
    Chapter 6 - Working with Nextflow schema
  </a>
</p>
<p class="text-left">
  <a href="linting_modules/" class="btn btn-lg btn-success" target="_blank">
    Chapter 7 - Linting your modules
  </a>
</p>

:::note{title="Key nf-core tools commands"}

- `nf-core create <pipeline>` creates a pipeline from the nf-core template.
- `nf-core lint` lints the pipeline code for things that must be completed.
- `nf-core modules list local` lists modules currently installed into your pipeline.
- `nf-core modules list remote` lists modules available to install into your pipeline.
- `nf-core modules install <tool/subtool>` installs the tool module into your pipeline.
- `nf-core modules create` creates a module locally to add custom code into your pipeline.
- `nf-core modules lint --all` lints your module code for things that must be completed.
- `nf-core schema build` opens an interface to allow you to describe your pipeline parameters and set default values, and which values are valid.

:::
