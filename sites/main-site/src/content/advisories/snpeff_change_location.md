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
      - "1.0.0"
      - "1.1.0"
      - "1.1.1"
      - "1.2.0"
      - "1.2.1"
  - name: sarek
    versions:
      - "3.0"
      - "3.0.1"
      - "3.0.2"
      - "3.1"
      - "3.1.1"
      - "3.1.2"
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

In August 2025, snpEff moved the hosting location of its caches.
snpEff versions older than 5.3 still point to the old location, so `snpEff download` fails with an error such as:

```console
FATAL ERROR: Failed to download database from [https://snpeff.blob.core.windows.net/databases/...
```

This affects every pipeline version that downloads the snpEff cache with an older snpEff, including sarek before 3.6.0 and rnavar before 1.2.2.

Newer snpEff versions download from the new location, but cache version 105 is no longer available there.

# Resolution

Update to sarek 3.6.0 or later, or rnavar 1.2.2 or later.
These versions use snpEff 5.3 and default to snpEff database version 99 or earlier.

If you need cache version 105, use a copy you downloaded before the change, and pass it with `--snpeff_cache`.
