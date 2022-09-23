---
title: Workflow specificity
subtitle: There should only be a single pipeline per data / analysis type.
menu:
  main:
    weight: 20
---

The nf-core community was founded to allow different groups to collaborate on pipelines instead of reinventing the same workflows in each institute.
As such, different pipelines should not overlap one another too much: there should only be a single pipeline for a given data + analysis type.
However, if the _purpose_ of the pipeline tasks and results are different, then this should be a separate pipeline.

If you would like to use a different set of tools to do a comparable analysis, then this should be added to the existing pipeline instead of creating something new.
