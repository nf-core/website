---
title: Software licences
subtitle: Automatically print software licences for DSL1 pipelines.
weight: 50
---

:::warning
This command does not currently work for DSL2 pipelines.
See [tools#905](https://github.com/nf-core/tools/issues/905)
and [tools#1712](https://github.com/nf-core/tools/issues/1712) for updates.
:::

Sometimes it's useful to see the software licences of the tools used in a pipeline.
You can use the `licences` subcommand to fetch and print the software licence from each conda / PyPI package used in an nf-core pipeline.

<!-- RICH-CODEX
timeout: 10
-->

![`nf-core licences deepvariant`](/images/tools/nf-core-licences.svg)
