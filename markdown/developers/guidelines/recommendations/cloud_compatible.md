---
title: Cloud compatible
subtitle: Pipelines should be tested on cloud computing environments.
menu:
  main:
    weight: 230
---

Pipelines should have explicit support for running in cloud environments.

Specifically, they should use the automated cloud-compute tests that are available for free through nf-core.
The [nf-core template](/tools#creating-a-new-workflow) comes with all required code to support this setup.

Where possible, hope that nf-core pipelines can also have benchmarks for pricing and performance when running on cloud environments.
