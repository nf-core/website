---
title: Trimming parameters ignored in RNA-seq pipeline 3.10
subtitle: The clip parameters for read trimming remain without affect for that pipeline version
category: pipelines
type: known_regression
severity: low
publishedDate: "2023-01-13"
reporter:
  - MatthiasZepper
reviewer:
  - MatthiasZepper
pipelines:
  - name: rnaseq
    versions:
      - "3.10.0"
      - "3.10.1"
modules:
subworkflows:
configuration:
nextflowVersions:
nextflowExecutors:
softwareDependencies:
references:
---

# Issue

Due to an oversight, the convenience parameters `--clip_r1`, `--clip_r2`, `--three_prime_clip_r1`, `--three_prime_clip_r2` that instruct Trim Galore to trim the respective number of bases from the reads are not used, even when specified by the user.
In rare cases, this may lead to slightly lower mapping rates.

# Resolution

This regression was not fixed.
The parameters were removed from the pipeline in version 3.11 to avoid confusion.
Since then, tool-specific trimming and other arguments can be provided via the `--extra_trimgalore_args` respectively the `--extra_fastp_args`.
