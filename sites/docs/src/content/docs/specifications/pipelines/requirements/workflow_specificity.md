---
title: Workflow specificity
subtitle: There should only be a single pipeline per data / analysis type.
menu:
  main:
    weight: 20
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

Different pipelines MUST NOT overlap one another: there MUST be only a single pipeline for a given data and analysis type.

If the purpose of the pipeline tasks and results are different, then this MUST be a separate pipeline.

Different sets of tools to perform comparable analysis SHOULD be added to the existing pipeline instead of creating a new pipeline.

:::info{title="Rationale" collapse}
The nf-core community was founded to allow different groups to collaborate on pipelines instead of reinventing the same workflows in each institute.
Maintaining a single pipeline per analysis type reduces duplication, consolidates community effort, and ensures better long-term maintenance and support.
:::

If in doubt, you can always ask on the [`#help`](https://nfcore.slack.com/channels/help) nf-core Slack channel, which is designed for exactly these kinds of discussions.
