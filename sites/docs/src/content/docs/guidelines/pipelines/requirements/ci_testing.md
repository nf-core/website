---
title: Continuous integration testing
subtitle: Pipelines must run CI tests.
shortTitle: CI testing
menu:
  main:
    weight: 80
---

Pipelines must have automated continuous integration testing, running using GitHub Actions.
There must also be CI tests using AWS (`test` and `test_full`).

There must be a config `profile` called `test` that should be as comprehensive as possible - that is, it should run as much of the pipeline as possible.
It should use as tiny test data set as possible (even if the output that it creates is meaningless).

## Increasing disk space for CI runners

nf-core uses [RunsOn](https://runs-on.com/) for self-hosted AWS runners. The default runner disk size may be insufficient for pipelines that generate large intermediate files or pull many container images. If your CI tests fail due to disk space, you can increase the volume size using the `volume` RunsOn job label.

```yaml title=".github/workflows/ci.yml"
jobs:
  test:
    runs-on:
      - runs-on=${{ github.run_id }}
      - runner=4cpu-linux-x64
      - volume=80gb
```

See the [RunsOn volume documentation](https://runs-on.com/configuration/job-labels/#volume) for full details.
