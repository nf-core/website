---
title: "nf-core/seqsubmit improvements"
category: pipelines
slack: https://nfcore.slack.com/channels/seqsubmit
location: Barcelona
image: "/assets/images/events/2026/hackathon-barcelona/seqsubmit_schema.png"
image_alt: "nf-core/seqsubmit workflow schema"
leaders:
  ochkalova:
    name: Sonya Ochkalova
    slack: https://nfcore.slack.com/team/U08AHEKTQQ4
  Vkale1:
    name: Varsha Kale
    slack: https://nfcore.slack.com/team/U0B6VM1PZV3
---

## Goal

[`nf-core/seqsubmit`](https://github.com/nf-core/seqsubmit) is a pipeline that automates submission of sequencing data to [ENA](https://www.ebi.ac.uk/ena/browser/home). We started this pipeline in October 2025 at the nf-core hackathon in Barcelona, and it recently got its first release (v1.0.0). At this hackathon we want to make it easier to test and contribute to, and extend it to cover cases it can't handle yet, such as co-assemblies, eukaryotic bins/MAGs and long read data.

## Description

The pipeline currently supports the following submission modes, each routed to a dedicated workflow:

- reads — raw sequencing reads submission via the READSUBMIT workflow (pink)
- metagenomic_assemblies — assembly submission via the ASSEMBLYSUBMIT workflow (green)
- mags — metagenome-assembled genomes (MAGs) submission via the GENOMESUBMIT workflow (blue)
- bins — bins submission via the GENOMESUBMIT workflow (blue)

Behind the scenes, the pipeline validates your input samplesheet, computes any missing statistics (coverage via CoverM, completeness/contamination via CheckM2, taxonomy via CAT_pack/BAT for MAGs and bins), builds the ENA-compatible manifest files each submission type needs, and submits everything to ENA using their own [webin-cli](https://github.com/enasequence/webin-cli) tool. See the [pipeline README](https://github.com/nf-core/seqsubmit#readme) for more on what it does and how it's structured.

The tasks below are a mix of things that make the pipeline friendlier to contribute to (like fixing the test suite so anyone can run it with their own ENA Webin account) and things that extend what it can submit (co-assemblies, eukaryotic bins/MAGs, long reads).

## Tasks

We probably won't get through all of these during the hackathon. We'll pick from the list below depending on how many people join and their experience level.

:::info{title="Great first issues"}
Several tasks below are Easy or Medium and can be done without prior nf-core experience, they're a good place to start if this is your first contribution to the nf-core pipeline.
:::

- Make the test suite runnable without maintainer Webin credentials ([#95](https://github.com/nf-core/seqsubmit/issues/95)) (Difficulty: Easy)
- Add tests for --is_private flag for different kinds of data ([#80](https://github.com/nf-core/seqsubmit/issues/80)) (Difficulty: Easy)
- Add `reads` mode output files to multiqc report ([#70](https://github.com/nf-core/seqsubmit/issues/70)) (Difficulty: Medium)
- Support submission of metagenomic assemblies generated from co-assemblies ([#61](https://github.com/nf-core/seqsubmit/issues/61)) (Difficulty: Medium)
- Support submission of MAGs/bins generated from co-assemblies ([#96](https://github.com/nf-core/seqsubmit/issues/96)) (Difficulty: Hard)
- Add eukaryotic MAGs/bins support ([#40](https://github.com/nf-core/seqsubmit/issues/40)) (Difficulty: Hard)
- Pipeline should support long-read submissions end-to-end ([#98](https://github.com/nf-core/seqsubmit/issues/98)) (Difficulty: Hard)

Each linked issue has background on why it matters, the terms you'll need, and a suggested step-by-step approach, so pick whichever one matches your interest and experience level.

## Making PRs

All updates should be made against the `dev` branch of the [nf-core/seqsubmit](https://github.com/nf-core/seqsubmit) repository.

Join our Slack channel [`#seqsubmit`](https://nfcore.slack.com/channels/seqsubmit) to get help or more information.
