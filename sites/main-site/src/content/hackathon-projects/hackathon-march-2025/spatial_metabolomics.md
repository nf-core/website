---
title: Bringing spatial metabolomics to nf-core
category: pipelines
intro_video: ""
slack: https://nfcore.slack.com/archives/C08JKPV4ZNH
image: "/assets/images/events/2025/hackathon-march/Spat_metabo_meme.jpg"
image_alt: "'Excited Jonah Hill' meme, where Jonah hill is screaming from excitment when he heard that spatial metabolomics workflows can be automated with nextflow"
leaders:
  Bisho2122:
    name: Bishoy Wadie
    slack: "https://nfcore.slack.com/team/U06JZ8F8LQG"
---

## Introduction

[METASPACE](https://metaspace2020.org/) is the cloud computing engine for annotating metabolites, lipids, and glycans in imaging mass spectrometry data, and hosts a knowledgebase of over 13.000 public annotated datasets from over 200 labs.

Most METASPACE users access the platform through the web app, which offers core functionalities for uploading, processing, and visualizing results.

However, bioinformaticians looking to integrate METASPACE into their workflows must develop custom scripts to interact with the METASPACE Python client and submit datasets.

Developing nf-core modules will streamline this process, enabling scalable and reproducible integration of METASPACE into existing bioinformatics pipelines.

This will also support METASPACE in-the-loop applications, particularly in spatial multi-omics integration workflows.

## Goal

1. Create the following nf-core modules using [Cardinal](https://www.bioconductor.org/packages/release/bioc/html/Cardinal.html):
   1. Data import (creating MSImagingExperiment object).
   2. Preprocessing (Primarily peak picking) to generate a centroided dataset -> input for METASPACE.

2. Create an API module to submit datasets to METASPACE using [METASPACE python client](https://metaspace2020.readthedocs.io/en/latest/)

3. Create a module to convert METASPACE datasets to AnnData and SpatialData objects using [METASPACE Converter](https://github.com/metaspace2020/metaspace-converter)

4. Create a module to perform enrichment analysis using [S2IsoMEr](https://github.com/alexandrovteam/S2IsoMEr) from AnnData object. (I'll personally work on this module)

5. If time permits, combine the above modules in one pipeline

We welcome contributors of all experience levels. All of the above goals are suitable for beginners.
