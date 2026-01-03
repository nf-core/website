---
title: Use Conda
subtitle: Package software using bioconda and biocontainers.
menu:
  main:
    weight: 200
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

All pipeline software SHOULD be packaged using either [Bioconda](https://bioconda.github.io/) or [conda-forge](https://conda-forge.org/).

:::info{title="Rationale" collapse}
Using conda for software packaging is strongly preferred because:

- Many pipeline users may be unable to use Docker / Singularity and run using Conda
- Conda packages are automatically packaged into Docker / Singularity images using [biocontainers](https://biocontainers.pro/) (Bioconda) or [Seqera Containers](https://seqera.io/containers/) (bioconda and conda-forge).
- We have automation to detect new bioconda / conda-forge software releases, which can trigger automated nf-core module updates

Note that all Bioconductor R packages are available via Bioconda.
Seqera Containers makes it easy to bundle multiple conda packages together into a single container, as well as providing images for both `aarch64` and `arm64` CPU architectures.
:::

Anyone can package software tools into Bioconda / conda-forge - you do not need to be the original author of the software.
Many nf-core developers have packaged software for conda, which benefits the wider bioinformatics community.

To get help with packaging on conda, please see the [`#bioconda` slack channel](https://nfcore.slack.com/archives/CM46YC6BZ).
