---
title: Workflow name
subtitle: Names should be lower case and without punctuation.
menu:
  main:
    weight: 35
---

All nf-core pipelines must have consistent naming that adhere to the following rules:

- Lower case only
- No special characters (a-z only)
- Descriptive and intuitive
- Must not clash with different data or analysis types

The first two points are to to maximise compatibility with other platforms,
so that the pipeline assets are named consistently across systems.

We strongly recommend that names are self-descriptive towards the data or analysis type the pipeline will be using or performing.
Users of the pipeline should be able to get a good idea of what the pipeline does from just the name.

Names must be approved by the nf-core community [(nf-core/proposals)](https://github.com/nf-core/proposals/issues).

Please refer to your pipeline with the `nf-core` prefix in written materials, for example `nf-core/yourpipeline`.

:::note
In _very_ rare occasions, we may be willing to bend this requirement. This has happened in the past when
a pipeline predates nf-core or comes with a substantial existing community.
This has happened for only a small handful of pipelines out of the >100 that we have in the community.
:::
