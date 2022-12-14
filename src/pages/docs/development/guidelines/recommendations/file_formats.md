---
title: File formats
subtitle: Use community accepted modern file formats.
menu:
  main:
    weight: 210
---

Pipelines should work with best practice modern file formats, as accepted by the community.

Where possible, genomics pipelines should generate `CRAM` alignment files by default, but have a `--bam` option to generate `BAM` outputs if required by the user.
