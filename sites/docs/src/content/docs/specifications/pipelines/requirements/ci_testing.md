---
title: Continuous integration testing
subtitle: Pipelines must run CI tests.
shortTitle: CI testing
menu:
  main:
    weight: 80
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

Pipelines MUST have automated continuous integration testing running using GitHub Actions.
There MUST also be CI tests using AWS (`test` and `test_full`).

There MUST be a config `profile` called `test` that SHOULD be as comprehensive as possible - that is, it SHOULD run as much of the pipeline as possible.
It SHOULD use as tiny a test data set as possible (even if the output that it creates is meaningless).
