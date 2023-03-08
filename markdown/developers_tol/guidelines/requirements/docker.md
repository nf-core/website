---
title: Docker Support
subtitle: Software must be bundled using Docker and versioned.
menu:
  main:
    weight: 70
---

Pipelines must have all software bundled using [Docker](https://www.docker.com/) - that is, it must be possible to run the pipeline with `-profile docker` and have all software requirements satisfied.

Tools should use docker images from [biocontainers](https://biocontainers.pro) where possible, as using Bioconda / Biocontainers gives support for conda + docker + singularity.

All containers must have specific, stable versions pinned.
These should preferably be named after a software release, but it can also be by commit or some other identifier.
Software versions must be static and stable. Labels such as `latest`, `dev`, `master` and so on are not reproducible over time and so not allowed.
