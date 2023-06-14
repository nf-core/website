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

If in doubt, you can always ask on the [`#new-pipelines`](https://nfcore.slack.com/archives/CE6SDEDAA) nf-core Slack channel, which is designed for exactly these kinds of discussions.
