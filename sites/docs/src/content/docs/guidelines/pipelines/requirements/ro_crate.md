---
title: Research Object Crate
subtitle: Pipelines must come with their own Research Object (RO) Crate
shortTitle: RO Crate
---

All nf-core pipelines MUST come with their own Research Object (RO) Crate.
RO Crate is an open-source standard which we use to describe our pipelines and their components in a structured way and helps with automated provenance tracking.

:::note
The RO Crate generated with nf-core/tools describes the pipeline as a whole, not pipeline runs.
For that kind of provenance, use the new [nf-prov plugin](https://github.com/nextflow-io/nf-prov), which is currently in development.
:::

More information can be found [here](https://www.researchobject.org/ro-crate/).
