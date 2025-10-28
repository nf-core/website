---
title: Docker Support
subtitle: Software must be bundled using Docker and versioned.
menu:
  main:
    weight: 70
---

Pipelines must have all software bundled using [Docker](https://www.docker.com/) - that is, it must be possible to run the pipeline with `-profile docker` and have all software requirements satisfied.

All containers must have specific, stable versions pinned.
Software versions must be static and stable. Labels such as `latest`, `dev`, `master` and so on are not reproducible over time and so not allowed.

:::note
Tools should be packaged with Bioconda / Conda-forge if possible.
Docker and Singularity images can then be automatically created.
See the ["Use Conda" recommendation](../recommendations/bioconda).

If packaging with Conda is impossible and a custom docker image _must_ be used,
please contact the nf-core `@core-team`.
See the ["Custom Docker images" recomendation](../recommendations/custom_containers)
:::
