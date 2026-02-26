---
title: nf-core/scrnaseq
category: pipelines
slack: https://nfcore.slack.com/archives/C06PKQ25BJT
location: Helmholtz Munich, Online
image: ""
image_alt: ""
leaders:
  Simone:
    name: Simone Röh
    slack: https://nfcore.slack.com/team/U07LWDWHU01
  Alena:
    name: Alena Boos
    slack: https://nfcore.slack.com/team/U07U6BW9W5U
---

Contribute to [nf-core/scrnaseq](https://nf-co.re/scrnaseq/dev/) pipeline.

## Goal

We aim to extend the nf-core/scrnaseq pipeline with functionality for RNA velocity calculation, sample aggregation, and alternative chemistries.

## Tasks

1. Integrate the existing `velocyto` module for RNA velocity calculation from Cell Ranger .bam outputs.
2. Create `cellranger_aggr` and `cellrangermulti_aggr` modules for Cell Ranger aggregation workflows.
3. Assess the feasibility of integrating the BD Rhapsody processing pipeline and, if viable, create a module for it. [Intro](https://www.bdbiosciences.com/en-us/products/software/rhapsody-sequence-analysis-pipeline) · [Manual](https://bd-rhapsody-bioinfo-docs.genomics.bd.com/top_introduction.html)
4. Assess the feasibility of integrating cell type annotation options and, if viable, create the necessary module or subworkflow. (matching 10X cellranger annotate)


:::note
Contributors of all experience levels are welcome; familiarity with Nextflow is helpful.
:::