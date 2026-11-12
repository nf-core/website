---
title: Merged families use the wrong sample's seed MSAs in proteinfamilies pipeline 2.0.0 - 2.4.0
subtitle: MERGE_FAMILIES:MERGE_SEEDS paired every sample with the first sample's seed alignments
category: [pipelines]
type: known_regression
severity: high
publishedDate: "2026-08-11"
reporter:
  - vagkaratzas
pipelines:
  - name: proteinfamilies
    versions: ["2.0.0", "2.1.0", "2.2.0", "2.3.0", "2.4.0"]
modules:
  - MERGE_SEEDS
subworkflows:
  - MERGE_FAMILIES
configuration:
nextflowVersions:
nextflowExecutors:
softwareDependencies:
references:
---

# Issue

The `MERGE_FAMILIES` subworkflow passed the seed MSA channel to `MERGE_SEEDS` as `seed_msa.first()`.
This converts the channel to a value channel holding only the **first** emission, so every pooled group
of similar families, from every sample, was merged against the first sample's seed alignments.

For a run with more than one input sample and family merging enabled (the default,
`--skip_family_merging false`), all merged families of the second and later samples were built from
sequences belonging to the first sample. Downstream `GENERATE_FAMILIES` outputs derived from those
merged seeds (MSAs, HMMs, recruited sequences and final family assignments) are therefore incorrect
for the affected samples.

Runs with a single input sample, and runs with `--skip_family_merging true`, are not affected.

# Resolution

`MERGE_FAMILIES` now combines the pooled components with the seed MSA collection by sample id, so each
pooled family is paired with its own sample's seeds while still allowing multiple groups per sample.
The fix is available in version 2.5.0 onwards ([#182](https://github.com/nf-core/proteinfamilies/pull/182)).

Users who ran an affected version on a multi-sample samplesheet without `--skip_family_merging` should
re-run with version 2.5.0 or later.
