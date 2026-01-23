---
title: Workflow name
subtitle: Names should be lower case and without punctuation.
menu:
  main:
    weight: 35
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

All nf-core pipelines MUST have consistent naming that adheres to the following rules:

- Lower case only
- No special characters (a-z only)
- Descriptive and intuitive
- MUST NOT clash with different data or analysis types

:::info{title="Rationale" collapse}
The first two points maximise compatibility with other platforms, so that the pipeline assets are named consistently across systems.
:::

Names SHOULD be self-descriptive towards the data or analysis type the pipeline will be using or performing.
Users of the pipeline SHOULD be able to get a good idea of what the pipeline does from just the name.

Names MUST be approved by the nf-core community [(nf-core/proposals)](https://github.com/nf-core/proposals/issues).

Pipelines SHOULD be referred to with the `nf-core` prefix in written materials, for example `nf-core/yourpipeline`.

:::note
In _very_ rare occasions, we may be willing to bend this requirement.
This has happened in the past when
a pipeline predates nf-core or comes with a substantial existing community.
This has happened for only a small handful of pipelines out of the >100 that we have in the community.
:::
