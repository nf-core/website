---
title: Cloud compatible
subtitle: Pipelines should be tested on cloud computing environments.
menu:
  main:
    weight: 230
---

Pipelines should have explicit support for running in cloud environments.

Specifically, they should use the automated cloud-compute tests that are available for free through nf-core.
The pipelines created with [nf-core template](https://nf-co.re/tools/#creating-a-new-pipeline) comes with all required code to support this setup.

Where possible, we hope that developers can create benchmarks for nf-core pipelines describing pricing and performance when running on cloud environments.
