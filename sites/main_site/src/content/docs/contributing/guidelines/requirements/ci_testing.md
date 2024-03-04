---
title: Continuous integration testing
subtitle: Pipelines must run CI tests.
menu:
  main:
    weight: 80
---

Pipelines must have automated continuous integration testing, running using GitHub Actions.
There must also be CI tests using AWS (`test` and `test_full`).

There must be a config `profile` called `test` that should be as comprehensive as possible - that is, it should run as much of the pipeline as possible.
It should use as tiny test data set as possible (even if the output that it creates is meaningless).
