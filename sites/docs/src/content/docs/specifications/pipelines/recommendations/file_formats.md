---
title: File formats
subtitle: Use community accepted modern file formats.
menu:
  main:
    weight: 210
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

Pipelines SHOULD work with best practice modern file formats, as accepted by the community.

Where possible, genomics pipelines SHOULD generate `CRAM` alignment files by default, but SHOULD have a `--bam` option to generate `BAM` outputs if required by the user.
