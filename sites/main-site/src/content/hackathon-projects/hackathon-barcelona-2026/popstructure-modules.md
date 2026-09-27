---
title: Population-genetics modules (TreeMix and PLINK)
category: components
slack: https://nfcore.slack.com/archives/C0C34GCCPN1
location: Barcelona
image: "/assets/images/events/2026/hackathon-barcelona/popstructure-modules-fry.png"
image_alt: "Meme of Fry from Futurama squinting suspiciously, captioned: Not sure if it's a bottleneck or just population structure. A strip below shows a small population tree with a dashed migration edge and the text: Let's find out with TreeMix + PLINK modules, nf-core hackathon Barcelona 2026."
leaders:
  malghuraybi:
    name: Mashael Alghuraybi
    slack: https://nfcore.slack.com/team/U0A9XPGSMSM
  NouraMalharbi:
    name: Noura Alharbi
    slack: https://nfcore.slack.com/team/U0BQUD29125
  Othmanaljurayyad:
    name: Othman Aljurayyad
    slack: https://nextflow.slack.com/team/U0BHKD48AFP
---

## Goal

Add missing building blocks for population-structure and population-history analysis to nf-core/modules: a TreeMix module, conversion of PLINK allele frequencies into TreeMix input, and a PLINK per-population allele-frequency module. Each should be merged, or at least in review, by the end of the hackathon, with `nf-test` coverage and a complete `meta.yml`.

## Description

nf-core/modules already covers several steps of genotype-based population-structure workflows, including PLINK and PLINK 2 modules for filtering, extraction, LD pruning and PCA (`plink/extract`, `plink/indep`, `plink2/filter`, `plink2/pca`, among others), as well as an `admixture` module.

What is currently missing are the components needed to extend these analyses toward population relationships and migration history with [TreeMix](https://bitbucket.org/nygcresearch/treemix/wiki/Home).

These modules are intended as reusable building blocks for population-genetics workflows using SNP-array and PLINK-compatible data. They also form part of the proposed `nf-core/popstructure` workflow ([proposal #170](https://github.com/nf-core/proposals/issues/170)), which overlaps in scope with the existing `consepopgen` pipeline.

The modules we plan to build:

| Module | What it does |
| --- | --- |
| PLINK allele frequencies | Calculates per-population allele frequencies from PLINK data (`--freq` with population assignments) |
| TreeMix input conversion | Converts per-population allele-frequency output into the format required by TreeMix |
| TreeMix | Fits a population tree with user-defined migration edges, root population, seed and block size |

Two design points we want to get right from the start:

- **Scientific choices stay explicit.** The number of migration edges, root population, block size and random seed are exposed as inputs rather than hidden defaults, so analytical choices are recorded and reproducible.
- **Species independence.** Non-standard chromosome sets (for example `--chr-set`) and population labels should work across organisms rather than assuming human chromosome conventions.

We checked nf-core/modules for existing work: there is currently no TreeMix module, and no PLINK module for per-population allele frequencies (`--freq`), runs of homozygosity (`--homozyg`) or genetic distances (`--distance`), with no open issue or PR for these components at the time of proposal.

## Tasks

Good project for anyone who wants to learn how nf-core modules are written, and for population geneticists who know these tools and want to see them wrapped properly.

:::info{title="Great first issues"}
This project is newcomer-friendly:
writing a small module and its `nf-test` is one of the best ways to learn nf-core module standards.
You will need a working Nextflow, `nf-core` tools and a container engine set up on your machine.
:::

Specific tasks:

- Open a module request issue for each missing module and link it to the project
- Write the PLINK per-population allele-frequency module
- Write the TreeMix input-conversion module
- Write the TreeMix module, including outputs for the tree, edges, residuals and likelihood files
- Add `nf-test`s with a small test dataset and snapshots that are stable across runs; TreeMix is stochastic, so a fixed seed will be used for testing
- Fill in `meta.yml`, `environment.yml` and topic channels for versions, then run `nf-core modules lint`
- Contribute a small population-genetics test dataset to nf-core/test-datasets if a suitable one does not already exist

### Stretch goals

- Add PLINK genetic-distance and runs-of-homozygosity modules
- Sketch a reusable population-structure subworkflow chaining the new modules, as a first step toward the proposed `nf-core/popstructure` pipeline

If you have data or a workflow that uses TreeMix, ADMIXTURE or PLINK, bring it. Real use cases will help us decide which parameters need to be exposed.

## Making PRs

Module PRs go against the `master` branch of [nf-core/modules](https://github.com/nf-core/modules).

Test data goes to [nf-core/test-datasets](https://github.com/nf-core/test-datasets).

Post in the project Slack channel to get paired with a reviewer, or find one of the leads at the table.
