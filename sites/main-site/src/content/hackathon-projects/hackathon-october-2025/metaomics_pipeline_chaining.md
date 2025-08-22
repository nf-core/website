---
title: meta-omics pipeline chaining
category: pipelines
intro_video: ""
slack: https://nfcore.slack.com/archives/C070369GP7T
image: "https://nf-co.re/logo/metaomics"
image_alt: "'nf-core/metaomics logo"
leaders:
  jfy133:
    name: James Fellows Yates
    slack: "https://nfcore.slack.com/team/UEM37TBAR"
  erikrikarddaniel:
    name: Daniel Lundin
    slack: "https://nfcore.slack.com/team/UR949F5LG
---

## Project Aim

Last year the [meta-omics special interest group](https://nf-co.re/special-interest-groups/meta-omics) started implementing a relatively simple system for pipeline chaining: a generic subworkflow for generating samplesheets for downstream pipelines, however we encountered a few issues that meant the project stalled.

We want to resurrect this project exploiting new Nextflow functionality that should make this much easier!

Anyone working on meta\*omics (metagenomics, metataxonomics, metatranscriptomics, metaproteomics, metabolomics) are welcome to join!

## Goals

1. Adapt previous samplesheet-generation implementation to use the new Nextflow publishing system
2. Implement samplesheet-generation across multiple nf-core/meta\*omics pipeliens
3. Add generic template to nf-core/tools
