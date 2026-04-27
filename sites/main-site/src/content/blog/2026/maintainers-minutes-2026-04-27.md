---
title: "Maintainers Minutes: March 2026"
subtitle: "Keeping you informed of the latest maintainers discussions"
pubDate: 2026-04-27T10:00:00+01:00
headerImage: "/assets/images/blog/maintainers-minutes-2024-07-26/maintainers-wide.png"
headerImageAlt: "Cartoon yellow rubber duck with nf-core logo badge on its body with the nf-core logo."
embedHeaderImage: false
authors:
  - "LouisLeNezet"
label:
  - "maintainers"
---

The 'Maintainers Minutes' aims to give further insight into the workings of the [nf-core maintainers team](/governance#maintainers) by providing brief summaries of the monthly team meetings.

## Overview

No, you didn't miss any maintainers' minutes; we missed publishing maintainers' minutes after November — apologies for the gap.
This post summarizes the main discussions from the four past maintainers' meetings held over the first months of 2026.

Many things happened during the first months of 2026, we had a hackathon, several lively debates on writing guidelines, deprecation, CI/infra, and contributor experience.

Major topics covered in this period:

- Writing guidelines
- Deprecation process
- March hackathon
- Miscellaneous things
- Bitesize talks

## Writing guidelines

We had lively discussions about relaxing or clarifying some writing rules to make modules and workflows easier to develop and maintain.

### Modules writing

The "one tuple to rule them all" school of thought (having only one input tuple) is gaining traction (big win for @Maxime !).
It simplifies parallelization and might be supported by Nextflow soon.
We will not actively rewrite all module inputs now, since adoption of record types will require further changes anyway.

Maintainers remain split on allowing new fields like `task.ext.prefix2`.
These can be pragmatic solutions for modules that genuinely need two strings (e.g., BWAmem index naming), but risk inconsistent inputs if used without documentation.
In conclusion to the debates, similarly to `task.ext.args2`, `task.ext.prefix2` can now be accepted if no simpler alternative is possible.

Another topic that sparked a lot of discussion was adopting a paired reference+index input channel.
Currently, many modules take the FASTA and its index in separate channels, which is far from ideal.
We now recommend keeping the FASTA paired with its index (same channel) when the module actually needs the index.
Whether a module requires the index can be tricky to check; if not explicit in the documentation, omitting the index and seeing if the module generates it can be sufficient (e.g., `samtools` modules).
We will not enforce a global channel standard across all modules — modules that do not need auxiliary files can be adapted in pipelines using a simple `.map{}` operator.

### Workflows writing

Emitting outputs from nested workflows is error-prone and hard to maintain because each level requires passing through identical outputs. This becomes cumbersome in pipelines with many emits (see rnaseq PR #1679 and the ongoing [Nextflow discussion](https://github.com/nextflow-io/nextflow/issues/6756)). We expect records/record-like features in Nextflow to help resolve this.

## Deprecation process

### Modules

The deprecation process prompted a long discussion after we deprecated `cat_cat` in favor of `find_concatenate`.
We want a lightweight, documented way to deprecate modules without leaving the modules repository cluttered with unused processes.

Agreed approach:

- Mark deprecated modules with a deprecation date and an expiry date.
- At the annual spring clean, delete modules that have been deprecated for more than one year.

We considered more advanced options (e.g., having nf-core tools infer deprecation from git history), but the infra team is currently capacity-constrained, so those remain future work.

There are also naming issues in the module repo: 39 process names don’t match the primary tool used in their script blocks (see [linked issue](https://github.com/nf-core/modules/issues/11073)). We should name modules to reflect the main tool and subtool (e.g., a module that primarily runs `AuthentiCT deam2cont` should be named `AUTHENTICT_DEAM2CONT` even if it calls `samtools view` beforehand).

### Subworkflows

Rike has done a large review of subworkflows and flagged many thin wrappers (one or two tools, little logic) as candidates for cleanup or deprecation (see [linked issue](https://github.com/nf-core/modules/issues/11075)).
Removing these will reduce maintenance overhead when updating widely used modules (looking at you, `samtools/*`).

## March Hackathon

The spring cleaning hackathon ran from 2–6 March and focused on tidying the modules repository ahead of the nf-core spring hackathon.

- ~100 issues and PRs were closed across the modules repo.
- Tags were added to improve triage and track work during and after the hackathon.
- Pipeline proposals were reviewed and maintainers were asked to add relevant issues to the hackathon project board.

Louis took over the post-hackathon cleanup for modules topics-migration opened PRs and issues (7 issues left).

## Miscellaneous

### Codespace / Ona

GitHub Codespaces (now Ona) launched an open-source plan.
Ana investigated organization-level options but hasn’t received a human response from GitHub yet.

### Meeting time

Our maintainer team now spans more time zones.
This is great for global participation but means some people can’t join current meeting times.
Options discussed include a monthly time shift or changing the meeting day — we’ll run a poll to decide.

### R nf-core utils package

R modules currently duplicate code to parse Nextflow inputs and produce log/version YAML files.
Louis created a small R package exposing two helper functions to handle these tasks.
It currently uses base R but could adopt argparse in the future.

### Bulk modules testing

Updating modules only to discover broken tests is painful.
We plan to run bulk tests on modules that haven't been run for a long time, produce a report of failures, and prioritize fixes.

### Singularity timeouts and new container system

Timeouts from Docker to Singularity conversions remain common.
Planned fixes will also produce conda lock files.
We may add a new command (e.g., `nf-core modules containers`) to better manage containers.

### nf-tests and strict syntax

Strict syntax support has been added to nf-test and released.

### Proprietary software

Some modules lack containers because the software is not open‑source.
In a few cases, vendors have provided licenses that allow nf‑core to run unit tests, but many remain untestable.

Example: `Dragen` requires dedicated hardware, so we can only run stub tests and cannot verify reported versions.

Therefore, all these modules should now have a prominent warning at the top, explaining they are untested or partially tested and that installation/use requires caution.

## Bitesize talk ideas

We have plenty of ideas but not enough volunteers! Here are some subjects we would like to cover:

- nf-metro
- nf-core <3 GPUs
- [Espanso plugin intended to help with reviewing PRs](https://hub.espanso.org/nf-core)
- [r-nfcore-utils package](https://github.com/nf-core/r-nf-core-utils)

## End

That's all folks!
We’ll aim to publish maintainers' minutes more regularly in the coming months.

\- :heart: from your #maintainers team!
