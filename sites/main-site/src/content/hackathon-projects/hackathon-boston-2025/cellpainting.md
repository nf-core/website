---
title: nf-core/cellpainting kickoff
category: pipelines
slack: https://nfcore.slack.com/archives/C08H9MWG06T
image: "/assets/images/events/2025/hackathon-boston/cell-painting.jpg"
image_alt: meme comparing traditional bioimage tools with OME tools for large scale processing
leaders:
  kenibrewer:
    name: Ken Brewer
    slack: https://nfcore.slack.com/team/U03MWA6LMV1
  maxulysse:
    name: Maxime Garcia
    slack: https://nfcore.slack.com/team/UE6D8290F
---

This project focuses on creating an nf-core pipeline for processing data for cell painting analysis using CellProfiler, cytotable, pycytominer, and other tools.

## Goals

- **Fill out the basic documentation for the planned pipeline in the nf-core/cellpainting repo.**
  - _For people who are familiar with cell painting and would like to focus on the documentation._
- **Add appropriate test data to nf-core/datasets for the pipeline.**
  - _For people who are familiar with cell painting and can identify data for use in unit tests._
- **Create a module(s) for CellProfiler processing steps of the workflow including:**
  - Z projection
  - QC
  - illumination correction
  - analysis
  - _For people who are familiar Nextflow and nf-core module development._

_We welcome contributors of all experience levels._
