---
title: BarcodeMiner
category: pipelines
# TODO: create a project-specific Slack channel and update this link
slack: https://nfcore.slack.com/channels/hackathon-mar-2027
location: in-person
image: "/assets/images/events/2027/hackathon-march/is-this-a-barcode-meme.jpg"
image_alt: "Is this a pigeon? meme, captioned 'is this a Barcode?'"
leaders:
  jorondo1:
    name: Jonathan Rondeau-Leclaire
    slack: https://nfcore.slack.com/team/U0BQNPDLRNZ
---

## Goal

Port an existing bash pipeline to Nextflow, write documentation. This will be my first time writing a Nextflow pipeline from the ground up. Anyone is welcome regardless of experience!

## Description

BarcodeMiner is a [pipeline in development](https://github.com/jflucier/envBarcodeMiner) that leverages Dicey ([Untergasser et al. 2026](https://doi.org/10.1093/nar/gkag445)) to mine the entire NIH Core_nt database for sequences that can be amplified by a given pair of primers. It allows building custom reference taxonomy database for amplicon-based sequencing experiments, but also identifying potential non-target organisms that a primer pair could amplify. It was first used in [Poirier et al, 2026](https://www.sciencedirect.com/science/article/pii/S0048969726002846) to build a taxonomic reference database for the trnL amplicon and is currently being used in a manuscript to profile arbuscular mycorrhizal fungal communities using the AML2-WANDA primer pair.

## Tasks

- Build a Dicey module
- Write subworkflows: database index builder, core algorithm, NCBI taxonomy cleaner
- Docs docs docs
- More on the way!
