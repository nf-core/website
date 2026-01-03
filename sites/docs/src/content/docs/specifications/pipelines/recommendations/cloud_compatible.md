---
title: Cloud compatible
subtitle: Pipelines should be tested on cloud computing environments.
menu:
  main:
    weight: 230
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

Pipelines SHOULD have explicit support for running in cloud environments.

Pipelines SHOULD use the automated cloud-compute tests that are available for free through nf-core.
The pipelines created with [nf-core template](/docs/nf-core-tools/pipelines/create) come with all required code to support this setup.

Developers SHOULD create benchmarks for nf-core pipelines describing pricing and performance when running on cloud environments where possible.
