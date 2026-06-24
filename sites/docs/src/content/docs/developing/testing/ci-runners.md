---
title: CI runners
subtitle: Configure RunsOn and fix CI disk space issues
shortTitle: CI runners
weight: 5
---

nf-core uses [RunsOn](https://runs-on.com/) for self-hosted AWS runners that run GitHub Actions.
For background on how nf-core CI is set up, see the blog post [State of the nf-core CI](https://nf-co.re/blog/2025/state-of-nf-core-CI).

## Increasing disk space for CI runners

The default runner disk size (30GB) may be insufficient for pipelines that generate large intermediate files or pull many container images.
If your CI tests fail due to disk space, you can increase the volume size using the `volume` RunsOn job label.

```yaml title=".github/workflows/nf-test.yml" {6}
jobs:
  test:
    runs-on:
      - runs-on=${{ github.run_id }}
      - runner=4cpu-linux-x64
      - volume=40gb
```

See the [RunsOn volume documentation](https://runs-on.com/configuration/job-labels/#volume) for full details.
