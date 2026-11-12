---
title: "CAPTn-seq: porting a Tn-seq analysis to Nextflow"
category: pipelines
slack: https://nfcore.slack.com/team/U0B46F0CQ74
location: Barcelona
image: "/assets/images/events/2026/hackathon-barcelona/captnseq-graphical-abstract.png"
image_alt: "Graphical representation of CAPTn-seq's broad steps, from read processing to gene essentiality calling and digenic interaction network generation"
leaders:
  simonjeanneau:
    name: Simon Jeanneau
    slack: https://nfcore.slack.com/team/U0B46F0CQ74
---

## Goal

Port a fully operational Snakemake pipeline to Nextflow, aiming for nf-core compliance, making it reproducible and portable.
This will allow other labs to process transposon insertion sequencing data, call gene essentiality and draw the genome-wide digenic interaction network (when used with a vast mutant collection) without rebuilding the analysis from scratch — a need the community genuinely has.

## Description

CAPTn-seq (Computational Analysis Pipeline for Tn-seq) is an end-to-end pipeline for transposon-insertion sequencing.
It covers two steps any Tn-seq study needs, for which no ready-made nf-core pipeline currently exists, namely 1) read processing and 2) per-gene essentiality calling.
On top of that it adds its core contribution: genome-wide digenic interaction calling, comparing relative essentiality across mutant backgrounds to identify synthetic gene-gene interactions directly from [High-Throughput Transposon Mutagenesis (HTTM)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10089323/) data, which is our variant of the Tn-seq methodology.

## Project Aim

The pipeline existed as a working but bespoke, single-lab Snakemake workflow, originally built for the Keio collection, using _E. coli_ strain BW25113 as its reference.
This project is porting it to Nextflow using nf-core conventions, specifically so it stops being "the pipeline that only runs on our cluster because one person knows all its assumptions" and becomes something any lab running Tn-seq or HTTM can clone and run.
The porting has already begun, but there is still custom code behind the original pipeline that needs to be containerized efficiently to reproduce the analysis cleanly, while adding all the extras required for nf-core compliance and hardening it into something genuinely reusable.

The interaction-calling stage is collection-agnostic, so it generalizes beyond the mutant collection it was first built on.

## Tasks

:::info{title="Disclaimer"}
The following list is only indicative.
I am actually quite new to Nextflow, so any help/guidance is more than welcome even if it is not an element of this list.
:::

- **(Medium)** Finish porting the remaining pipeline stages to Nextflow: sample-level reporting and the genome-wide digenic interaction calling (the pipeline's core contribution) are the two stages not yet in Nextflow
- **(Easy)** Swap hand-rolled local modules (cutadapt, BWA, samtools-based steps) for official nf-core/modules equivalents where they already exist upstream, which means increased maintainability
- **(Medium)** Validate every ported stage against the original Snakemake pipeline's output, not just "it runs without crashing". This is mechanical but requires judgment to tell real bugs apart from expected numeric noise
- **(Medium)** Reach nf-core compliance: pass `nf-core lint`, add nf-test coverage
- **(Medium)** Build a public test dataset and CI test profile: real sample data can't be shared, so there's currently no way to run this pipeline's tests without private data. This is needed for both CI and for other labs to trust the pipeline works before pointing it at their own data
- **(Easy/Medium)** Generalize site-specific assumptions into a documented, swappable config pattern so another lab can actually adopt this on their own infrastructure
- **(Easy)** Write real usage and output documentation: no extensive pipeline-internals knowledge required and a good stress test for the code readability
- **(Hard, stretch)** Pursue official nf-core pipeline submission once the above is solid enough to survive community review

Success: a working CAPTn-seq run that efficiently reproduces the original HTTM analysis on test data.

## Making PRs

PRs go against the [simonjeanneau/captnseq](https://github.com/simonjeanneau/captnseq) repository, the Nextflow draft of the pipeline.

There is also a [talk at Nextflow Summit 2026](https://summit.nextflow.io/2026/virtual/agenda/httm-bringing-genome-wide-digenic-interaction-screens-to-nf-core/) on October 13th about porting this pipeline to Nextflow.
