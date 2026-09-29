---
title: Metaboigniter updates and benchmarking
category: pipelines
slack: https://nfcore.slack.com/channels/metaboigniter
location: Barcelona
image: "/assets/images/events/2026/hackathon-barcelona/metaboigniter.jpg"
image_alt: "Photo showing 'project update' on Scrabble stones. Photo by Matilda Alloway (@matildaonthemove) on Unsplash."
leaders:
  enryh:
    name: Henry Webel
    slack: https://nfcore.slack.com/archives/C010AEBQ599
---

## Goal

Metaboigniter updates and evaluation of adding mzmine as a backend for computations.

## Description

[nf-core/metaboigniter](https://github.com/nf-core/metaboigniter) is nf-core's pipeline
for pre-processing mass-spectrometry metabolomics data. The last release,
[2.0.1](https://nf-co.re/metaboigniter/2.0.1/), is now over two years old. The pipeline
still seems to be used, but not maintained.

MzMine is a popular open-source alternative to the current OpenMS backend with a
possibility of batch processing. In principle users should be able to choose a backend
and also benchmark different backends against each other for their applications. Adding
common benchmark datasets is something that is still missing for metaboigniter as far
as I can tell.

## Tasks

We'll start by discussing the most urgent pipeline updates and see which expertise on
Metabolomics data processing we have in the room.

### Open Current issues:

- [sample-order mismatch between the ConsensusXML and quantification table
  (#113, Jan 2026)](https://github.com/nf-core/metaboigniter/issues/113),

- [add a step to ensure
  indexed-mzML files (#95, open since May 2024)](https://github.com/nf-core/metaboigniter/issues/95),

- [check on updates for PYOPENMS_MSMAPPING module (#103, Dec 2024)](https://github.com/nf-core/metaboigniter/issues/103).

- update to newer OpenMS releases

### Benchmarking MzMine on published benchmarking datasets:

- discuss pitfalls of benchmarking metabolomics data
- get an idea of good parameters for MzMine and OpenMS for benchmarking

### Docker image and template updates:

- OpenMS image is using OpenMS 3.0.0, where we by now are at version 3.5.0
- local models use python scripts from `bin` folder
- merge [3 PRs](https://github.com/nf-core/metaboigniter/pulls) with updates
