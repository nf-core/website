---
title: Adding removal of genome contaminants module to nf-core/sarek
category: pipelines
intro_video: ""
slack: https://nfcore.slack.com/archives/C05V9FRJYMV
image: "/assets/images/events/2025/hackathon-march/genome_contaminants.jpg"
image_alt: "Fry wondering if he'll win a Nobel Prize"
leaders:
  apsteinberg:
    name: Asher Preska Steinberg
    slack: "https://nfcore.slack.com/team/U07ABAY61SS"
---

This project will work on adding a module to filter out genomic contaminants to the [nf-core/sarek](https://nf-co.re/sarek/) pipeline. This would be analogous to what is currently implemented in the [nf-core/rnaseq](https://nf-co.re/rnaseq/3.18.0/) pipeline with [BBSplit](https://nf-co.re/rnaseq/3.18.0/docs/output/#bbsplit). We will focus on implementing [xengsort](https://gitlab.com/genomeinformatics/xengsort) to do this. The steps here would be:

- finishing off the xengsort module
- making a subworkflow to run xengsort
- integrating this subworkflow into sarek

All skill levels welcome.

## Goals

Introduce the filtering of genomic contaminants as a feature in the sarek pipeline.
