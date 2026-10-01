---
title: nf-core/mag improvements
category: pipelines
slack: https://nfcore.slack.com/channels/mag
location: Barcelona
leaders:
  jorondo1:
    name: Jonathan Rondeau-Leclaire
    slack: https://nfcore.slack.com/team/U0BQNPDLRNZ
---

## Goal

Improvements to the [nf-core/mag](https://github.com/nf-core/mag) pipeline.

## Description

The developers of nf-core/mag list several tools that could be added to this pipeline (see [nf-core/mag#1022](https://github.com/nf-core/mag/issues/1022)), as well as potential replacements — for example, QUAST is currently overkill and generates tons of useless files.
I'm but a humble recent contributor to that pipeline, and a novice at Nextflow, but I'd love to work on this with other motivated people.

Hackathon attendees could choose which of these to work on and tackle some of them.
There is also a discussion to be had regarding the feasibility and pertinence of an additional entry point (providing bins directly, for bin refinement), but this would require experienced Nextflow developers, which I cannot say that I am!

## Tasks

- Pick tools from [nf-core/mag#1022](https://github.com/nf-core/mag/issues/1022) to add to the pipeline
- Look at replacing tools that no longer pull their weight (e.g. QUAST)
- Discuss the feasibility of an additional entry point that takes bins directly, for bin refinement

:::note
Keep an eye on this page for more updates.
:::
