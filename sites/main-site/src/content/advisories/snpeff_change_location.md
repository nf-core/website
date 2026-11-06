---
title: snpEff cache location change in pipelines using snpEff
subtitle: All pipelines using snpEff are affected by the change of the default cache location in snpEff happening in August 2025.
category:
  - pipelines
  - configuration
type: incompatibility
severity: high
publishedDate: "2025-10-22"
reporter:
  - maxulysse
  -
reviewer:
  - maxulysse
  -
pipelines:
  - name: sarek
    versions:
      - 2.12.0
  - name: rnavar
    versions:
      - 4.2.0
modules:
  - name: snpeff_download
  - name: snpeff_snpeff
subworkflows:
  - name: vcf_annotation_snpeff
nextflowVersions:
  - ">=25.04.0 <25.04.3"
references:
  - title: snpeff issue regarding change of hosting location of cache and binaries
    description: Issues reporting the issue with snpeff download after the change of hosting location
    url: https://github.com/pcingola/SnpEff/issues/602
---

# Issue

https://github.com/pcingola/SnpEff/issues/602

snpeff download is broken with older versions of sarek/rnavar due to changes in location for the storage of snpeff cache.
snpeff download works with newer versions of sarek/rnavar with an updated of snpeff, but cache 105 is now not valid anymore (hopefully until a fix comes).

-> https://github.com/pcingola/SnpEff/issues/609
