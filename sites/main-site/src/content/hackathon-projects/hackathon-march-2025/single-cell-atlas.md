---
title: "Adapt Single Cell Expression Atlas analysis workflow to nf-core standards"
category: pipelines
intro_video: ""
slack: https://nfcore.slack.com/archives/C08J6LDLN0P
image: "/assets/images/events/2025/hackathon-march/single-cell-atlas.png"
image_alt: EMBL-EBI Single Cell Expression Atlas logo
leaders:
  pmb59:
    name: Pedro Madrigal
    slack: "https://nfcore.slack.com/team/U0447DB547L"
  irisdianauy:
    name: Iris Yu
    slack: "https://nfcore.slack.com/archives/D04P3TP4E1L"
---

The Single-cell Expression Atlas analysis workflow runs downstream analyses using Scanpy, leveraging the [scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts) package to run individual steps of the Scanpy workflow.

See [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc)
and [current Nextflow workflow](https://github.com/ebi-gene-expression-group/scxa-tertiary-workflow).

## Goal

Adapting the Single-Cell Expression Atlas Analysis Workflow to meet `nf-core` standards.
