---
title: nf-core/scalpel
category: pipelines
slack: https://nfcore.slack.com/archives/C0AST9JDVHC
image: "/assets/images/events/2026/hackathon-boston/SCALPEL_meme.png"
image_alt: "A photo of a surgeon bending over with a scalpel in the hand. Above the following dialogue is written: A: 'Relax Marcel, it's just a small cut with a scalpel, don't be nervous.' B: 'Doctor, my name isn't Marcel!' 'I know, I am Marcel.'"
leaders:
  Franzx7:
    name: Franz AKE
    slack: "https://nfcore.slack.com/team/U0888QNDQAJ"
---

This project focuses on improving [nf-core/scalpel](https://github.com/plasslab/nf-core-scalpel), a complex nf-core pipeline for large omics data with single-cell data, isoform-aware quantification, and downstream analysis. SCALPEL was previously developed as a standalone Nextflow pipeline and the associated work was published in [Nature Communications (2025)](https://www.nature.com/articles/s41467-025-61118-0). The current effort consists of reimplementing SCALPEL as a fully fledged nf-core pipeline so that it can be more easily deployed, maintained, and extended while integrating important updates into the workflow.

The Boston hackathon is a great opportunity to expand the pipeline scope, improve reliability, and make it easier for new contributors to get involved.

## Goal

By the end of the hackathon, we want to make clear progress towards a more robust and release-ready nf-core implementation of SCALPEL, with broader support for different single-cell input technologies and stronger downstream analysis capabilities.

## Planned tasks

- Integrate support for processing either 10x Chromium, Drop-seq, and potentially additional single-cell input formats.
- Work on nf-tests to improve coverage, reproducibility, and confidence in pipeline behavior.
- Set up sensible default configurations and container definitions for easier deployment.
- Integrate downstream single-cell analysis steps into the Nextflow execution.
- Add clearer error and warning messages for common or specific user-facing issues.
- Improve documentation and usage examples.

## Who should join?

This project is a good fit for people interested in:

- working on a complex Nextflow pipeline for large omics data processing, particularly single-cell data for isoform quantification
- developing pipelines to process different types of omics input files
- pipeline testing, reproducibility, and nf-test
- R or Nextflow-based downstream analysis development

## Recommended preparation

Participants will benefit from:

- Basic familiarity with Git and GitHub
- Some exposure to Nextflow or nf-core
- An interest in testing, debugging, or improving documentation
