---
title: Adapt SingleCell perturb-Seq CRISPR  pipelineto nf-core Standards and Implement Testing
category: pipelines
slack: https://nfcore.slack.com/archives/C0ARV562CNN
leaders:
  Lucas Ferreira da Silva:
    name: Lucas Ferreira da Silva
    slack: "https://nextflow.slack.com/team/U09F6344B5F"
---

![Alt text](https://github.com/pinellolab/CRISPR_Pipeline/blob/main/images/crispr_pipeline.png)

# Adapt CRISPR_Pipeline to nf-core Standards and Implement Testing

## Introduction

The [CRISPR_Pipeline](https://github.com/pinellolab/CRISPR_Pipeline) is a powerful tool for analyzing CRISPR screen data, but to maximize its utility, reproducibility, and integration with the broader bioinformatics community, it needs to align with [nf-core guidelines](https://nf-co.re/docs/guidelines).

Some of the modules within the current architecture are already compatible with nf-core standards. This project will focus on identifying those components, bridging the gaps with the remaining custom modules, and most importantly, establishing a robust, automated testing framework. To achieve meaningful testing without overloading CI/CD resources, we will implement a strategy to capture, downsample, and utilize public datasets—specifically targeting the Perturb-seq data from Gasperini et al. (2019) [Downloading the Pilot Dataset](#downloading-and-unpacking-a-small-fraction-of-the-gasperini-2019-pilot-dataset)
.

## Goals

1. **Audit and Align Modules:**
   - Review the existing `pinellolab/CRISPR_Pipeline` codebase to identify which Nextflow modules are already nf-core compatible.
   - Replace or refactor non-compliant components using nf-core `modules` and `subworkflows` where possible to ensure standardized genomic sequence modeling and pipeline execution.

2. **Dataset Curation for Testing:**
   - Retrieve the public Perturb-seq dataset from Gasperini et al. (2019).
   - Create a standardized script or Nextflow subworkflow to appropriately downsample this dataset so it is small enough for rapid CI/CD runs while retaining enough complexity to validate the pipeline's analytical steps.

3. **Implement CI/CD Testing Framework:**
   - Set up GitHub Actions workflows following nf-core templates.
   - Integrate the downsampled Gasperini dataset as the standard `test` profile.
   - Ensure that the test suite runs end-to-end, validating both single-module execution and the complete workflow.

4. **Documentation and Templating:**
   - Update the pipeline documentation (e.g., `usage.md`, `output.md`) to reflect the new nf-core structure and parameters.
   - Ensure all metadata and schema validations (`nextflow_schema.json`) are correctly implemented.

## Expected Outcomes

- A modernized, nf-core-compliant version of the CRISPR_Pipeline ready for either direct submission to nf-core or deployment as a highly standardized institutional tool.
- A reproducible, lightweight test profile using real-world downsampled Perturb-seq data.
- Successful, passing CI tests demonstrating that the refactored pipeline yields accurate and expected results.

## Getting Involved

If you are interested in Nextflow, CRISPR screens, or building robust bioinformatics infrastructure, we would love your help! Whether you want to focus on module refactoring, scripting the downsampling pipeline for the Gasperini dataset, or configuring GitHub Actions, there is plenty to tackle.

Drop into the `#boston-2026-crispr-pipeline` Slack channel to say hi and coordinate, or send a direct message to my [Nextflow Slack profile](https://nextflow.slack.com/team/U09F6344B5F)!
