---
title: Pass lint tests
subtitle: The pipeline must not have any failures in the `nf-core pipelines lint` tests
menu:
  main:
    weight: 130
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

The pipeline MUST NOT have any failures in the `nf-core pipelines lint` tests.

The `nf-core pipelines lint` command MUST be run by continuous integration tests on GitHub Actions and MUST be passing before release.

:::info{title="Rationale" collapse}
In order to automate and standardise the nf-core best practices, we have built a code linting tool.
These tests are run by the [nf-core/tools](https://github.com/nf-core/tools) package to ensure consistent quality and adherence to community standards across all pipelines.
:::

You can see the list of tests and how to pass them on the [error codes page](https://nf-co.re/tools).

In exceptional circumstances, certain tests MAY be ignored using `nf-core.yml`, if agreed upon by the community.
