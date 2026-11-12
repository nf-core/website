---
title: snpEff cache download fails in sarek and rnavar
subtitle: The snpEff cache moved to a new hosting location, so older pipeline versions can no longer download it.
category:
  - pipelines
  - modules
  - subworkflows
type: incompatibility
severity: high
publishedDate: "2026-11-12"
reporter:
  - maxulysse
reviewer:
  - mashehu
pipelines:
  - name: rnavar
    versions:
      - "1.1.0"
      - "1.1.1"
      - "1.2.0"
      - "1.2.1"
  - name: sarek
    versions:
      - "3.2.0"
      - "3.2.1"
      - "3.2.2"
      - "3.2.3"
      - "3.3.0"
      - "3.3.1"
      - "3.3.2"
      - "3.4.0"
      - "3.4.1"
      - "3.4.2"
      - "3.4.3"
      - "3.4.4"
      - "3.5.0"
      - "3.5.1"
modules:
  - snpeff/download
  - snpeff/snpeff
subworkflows:
  - vcf_annotate_snpeff
configuration:
nextflowVersions:
nextflowExecutors:
softwareDependencies:
references:
  - title: SnpEff issue 602
    description: GitHub issue about the problem with snpeff download after the change of hosting location
    url: https://github.com/pcingola/SnpEff/issues/602
  - title: SnpEff issue 609
    description: GitHub issue about snpEff cache version 105 no longer being available
    url: https://github.com/pcingola/SnpEff/issues/609
---

# Issue

In August 2025, snpEff moved the hosting location of its caches ([SnpEff#602](https://github.com/pcingola/SnpEff/issues/602)).
snpEff versions older than 5.3 still point to the old location, so downloading a cache with `--download_cache` fails with an error such as:

```console
FATAL ERROR: Failed to download database from [https://snpeff.blob.core.windows.net/databases/...
```

In addition, cache version 105 is no longer available at the new location ([SnpEff#609](https://github.com/pcingola/SnpEff/issues/609)).
The affected pipeline versions use version 105 by default for `GRCh38`, `R64-1-1` and `WBcel235`.

Affected versions, when run with `--download_cache`:

- sarek 3.2.0 to 3.5.1
- rnavar 1.1.0 to 1.2.1

Runs that use the default annotation cache or a cache you downloaded yourself are not affected.

# Resolution

Update to sarek 3.6.0 or later, or rnavar 1.2.2 or later.
These versions use snpEff 5.3 and default to cache version 99 instead of 105.

If you need cache version 105, do not use `--download_cache`.
The default annotation cache (`--snpeff_cache s3://annotation-cache/snpeff_cache/`) still provides it.
