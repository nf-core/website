---
title: Pass lint tests
subtitle: The pipeline must not have any failures in the `nf-core lint` tests
menu:
  main:
    weight: 130
---

In order to automate and standardise the nf-core best practices, we have built a code linting tool.
These tests are run by the [nf-core/tools](https://github.com/nf-core/tools) package.
The `nf-core lint` command must be run by continuous integration tests on GitHub Actions and must be passing before release.

You can see the list of tests and how to pass them on the [error codes page](https://nf-co.re/tools).

In some exceptional circumstances, it is ok to ignore certain tests using `nf-core.yml`, if agreed upon by the community.
