---
title: "The Mass Spectrometry Proteomics Special Interest Group"
subtitle: "Unifying efforts for mass spectrometry data analysis in nf-core"
pubDate: 2025-11-26T12:00:00+01:00
headerImage: https://images.unsplash.com/photo-1483428400520-675ef69a3bc4
headerImageAlt: "Photo by Hudson Hintze on Unsplash"
authors:
  - "DongzeHe"
tags:
  - proteomics
  - mass-spectrometry
  - special-interest-group
label:
  - "massspec-proteomics"
embedHeaderImage: false
---

We are excited to announce the formation of `massspec-proteomics`, the **Mass Spectrometry Proteomics** special interest group within nf-core!

## Introduction

As the use of Nextflow for proteomics data analysis continues to grow, so does the number of pipelines and modules being developed by the community. From established pipelines like `quantms` to newer initiatives like `mspepid` and `WOMBAT-P`, there is a wealth of expertise distributed across different teams. However, this growth also brings a challenge: how do we ensure we aren't reinventing the wheel?

The new MassSpec-Proteomics SIG has been created to serve as a central hub to coordinate these efforts. Our goal is to establish a cohesive, modular ecosystem in nf-core where components can be shared and pipelines can be interoperable.

## Our Vision

The group's primary focus is to establish guidelines and standards tailored for mass spectrometry data analysis. We are moving towards a vision where specific data types, such as **DDA** (Data-Dependent Acquisition), **DIA** (Data-Independent Acquisition), and **TMT** (Isobaric Labeling), are handled by dedicated, thin, and modular pipelines that can be easily chained together.

To achieve this, we are prioritizing:

1. **Shared Components:** Systematically refactoring and contributing high-quality modules and subworkflows (e.g., from `quantms` and for tools like FragPipe) to the official `nf-core/modules` repository.

2. **Standardization:** Establishing best practices for QC and benchmarking to ensure all proteomics pipelines meet the same high standards.

3. **Collaboration:** providing a forum for developers to align on roadmaps and technical solutions.

## Join Us!

Whether you are a developer working on a specific tool, a bioinformatician building pipelines, or a user with feedback on current workflows, we want to hear from you.

- **Slack:** Join the discussion in the `#massspec-proteomics` channel on the [nf-core Slack](https://nf-co.re/join).

- **Meetings:** We will be organizing regular virtual meetings to discuss technical challenges and roadmap planning. Keep an eye on the Slack channel for the schedule.

We look forward to building the future of nf-core proteomics with you!
