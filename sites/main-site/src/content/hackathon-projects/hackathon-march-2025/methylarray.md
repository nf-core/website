---
title: methylarray
category: pipelines
intro_video: ""
slack: https://nfcore.slack.com/archives/C06PZ8FURPW
image: "/assets/images/events/2025/hackathon-march/methylarray.png"
image_alt: "nf-core/methylarray"
leaders:
  ajandria:
    name: Adrian Janucik
    slack: "https://nfcore.slack.com/team/U04MRN5KYB1"
  ghada:
    name: Ghada Nouairia
    slack: "https://nfcore.slack.com/team/U033PMQJC6M"
  rafalcode:
    name: Ramon Fallon
    slack: "https://nfcore.slack.com/team/U02UPM988LR"
---

A pipeline for processing Illumina DNA Methylation array datasets: 450k, EPIC and EPICv2

## Goals

1. tidy up whatever has been written now and write custom containers for the whole pipeline (3 processes for now)
2. attach the Dockerfile to the pipeline itself (cellranger in scrnaseq cannot be run with conda so I think this pipeline could rely on containers only until we find a fix)
3. scripts rely on the directory access/path rather than the file access/paths (difficulty in fetching the directory itself from the test-data repo: perhaps add module for creating the CSV file with the input path as the scripts currently are expecting it and align with the requirements for proper -profile test setup.
4. tidy up the R scripts and parametrize hardcoded options
5. add a disclaimer that this will not run with conda
6. update docs
7. coordinate if the output of the pipeline generates same results using previously analyzed projects (align with @Ghada Nouairia)
8. remove FastQC module and remove MultiQC module or add custom QC plots with versions to the report
9. currently two tests are failing on GitHub actions related to pipeline logos, but these pass the nfcore tools repo.

## Motivation

Differential analysis of DNA methylation can come from sequencing reads, but also from Illumina methylation array which
are for detecting methylation intensity at certain pre-programmed sites. This pipeline is for the latter technique as there
already is a nf-core project for the former (methylseq).
