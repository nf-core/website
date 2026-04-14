---
title: sarek parabricks integration
category: pipelines
slack: https://nfcore.slack.com/archives/C0AS4MC7UG7
leaders:
  gburnett-nvidia:
    name: Gary Burnett
    slack: https://nfcore.slack.com/team/U05LY7TP5SL
---

This hackathon project focuses on updating
the [Sarek](https://nf-co.re/sarek/3.8.1/) pipeline with more Parabricks GPU modules as optional replacements for CPU-based steps. This will reduce the total runtime of Sarek for select workflows. There is a [GitHub Issue](https://github.com/nf-core/sarek/issues/1853) open on this topic which details previous work.

---

## Goal

Add more GPU modules to the Sarek pipeline.

Presently, Sarek utilizes the Parabricks [fq2bam module](https://nf-co.re/modules/parabricks_fq2bam/) as a GPU accelerated version of BWA-Mem. However, there are other Parabricks modules that can be added to Sarek:

- [Haplotypecaller](https://nf-co.re/modules/parabricks_haplotypecaller)
- [Deepvariant](https://nf-co.re/modules/parabricks_deepvariant)
- [Mutectcaller](https://nf-co.re/modules/parabricks_mutectcaller)

---

## What participants will do

Each contributor will:

1. Choose a module to work on with the team
2. Add it to the Sarek pipeline
3. Run tests and lint checks locally.
4. Open a Pull Request for review.

Each module will have its own branch and PR.

---

## Recommended preparation

Participants should ideally have:

- Basic familiarity with Git and GitHub (forking repositories, creating branches, and opening pull requests)
- Basic knowledge of Nextflow
- Familiarity with nf-core modules and pipelines

The following training material is recommended:

- [Adding modules to pipelines](https://nf-co.re/docs/tutorials/nf-core_components/adding_modules_to_pipelines)
