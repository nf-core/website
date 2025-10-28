---
title: extended short reads QC subworkflow
category: components
intro_video: ""
slack: https://nfcore.slack.com/archives/C070369GP7T
image: "/assets/images/events/2025/hackathon-barcelona/nfcore-metaomics_logo.png"
image_alt: "nf-core/metaomics logo"
leaders:
  vagkaratzas:
    name: Evangelos Karatzas
    slack: "https://nfcore.slack.com/team/U05LNHCFLCW"
---

## Project Aim

In the [meta-omics special interest group](https://nf-co.re/special-interest-groups/meta-omics) we started a discussion about having a standardized, general subworkflow for reads quality check (QC) in nf-core.
Currently, different nf-core pipelines (and not only) use their own implementations of essentially the same QC steps: stats, adapter clipping, paired-end merging, umi-detection deduplication etc.

For example, some modules are shared across multiple pipelines, and are used widely (e.g., fastqc -73 pipelines, fastp -22 pipelines) while others are used in fewer pipelines (e.g., trimgalore -15 pipelines, umitools_extract -10 pipelines,adapterremoval -2 pipelines, etc.).

The proposed project envisions a single Illumina short read preprocessing subworkflow that is installed and consistent across all nf-core pipelines using short read DNA data. However, at the same time the subworkflow will provide users the freedom to pick their preferred tool from a pool of nf-core QC modules or completely skip some of the QC steps.

It will also function as a template for an equivalent potential long reads QC subworkflow.

More information can be found on this [link](https://hackmd.io/@nf-core/S1NVXogilg).

Anyone working with Illumina short read data, and particularly within meta-omics (metagenomics, metataxonomics, metatranscriptomics, metaproteomics, metabolomics), are welcome to join!

## Goals

1. Push any missing nf-core modules (i.e., seqfu_check).
2. Create the required sub-subworkflows (i.e., shortread_adapter_removal, shortread_host_removal, etc.) and push to nf-core.
3. (minimal success) Create a minimal short reads QC subworkflow that does the following; (i) inital QC + stats, (ii) barcoding, (iii) adapter removal + run merge, (iv) complexity filtering, (v) deduplication, (vi) host removal, and (vii) final QC, with the minimum required tools (1 module per task).
4. (success) Create an extended version of the subworkflow with all available QC tools for short reads; multiple module options per task to choose from + skip parameters.
5. (bonus) Initialize a similar subworkflow but for long reads QC.
