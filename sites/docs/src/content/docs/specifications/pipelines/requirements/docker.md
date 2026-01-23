---
title: Docker Support
subtitle: Software must be bundled using Docker and versioned.
menu:
  main:
    weight: 70
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

Pipelines MUST have all software bundled using [Docker](https://www.docker.com/) - that is, it MUST be possible to run the pipeline with `-profile docker` and have all software requirements satisfied.

All containers MUST have specific, stable versions pinned.
Software versions MUST be static and stable.
Labels such as `latest`, `dev`, `master` and so on are not reproducible over time and MUST NOT be used.

:::note
Tools SHOULD be packaged with Bioconda / Conda-forge if possible.
Docker and Singularity images can then be automatically created.
See the ["Use Conda" recommendation](../recommendations/bioconda).

If packaging with Conda is impossible and a custom docker image MUST be used,
please contact the nf-core `@core-team`.
See the ["Custom Docker images" recomendation](../recommendations/custom_containers)
:::
