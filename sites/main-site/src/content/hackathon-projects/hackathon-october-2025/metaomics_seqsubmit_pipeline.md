---
title: "SeqSubmit: Automated Data Submission"
category: pipelines
intro_video: ""
slack: https://nfcore.slack.com/archives/C09F9HR0A9E
image: "/assets/images/events/2025/hackathon-barcelona/nfcore-metaomics_logo.png"
image_alt: "nf-core/metaomics logo"
leaders:
  KateSakharova:
    name: Ekaterina Sakharova
    slack: https://nfcore.slack.com/team/U04LFFXEW6N
  ochkalova:
    name: Sonya Ochkalova
    slack: https://nfcore.slack.com/team/U08AHEKTQQ4
---

## Project Aim

Sequencing experiments generate valuable data that should be shared with the scientific community, encouraging [FAIR principles](https://www.go-fair.org/fair-principles/). Yet submission through public repositories remains a major bottleneck. Researchers must navigate a complex multi-step process to decipher schemas, lengthy documentation and manually organise metadata according to strict specifications—all before a single read file reaches the repository.
This complexity means that **data submission is universally dreaded and avoided** depriving the community of potentially valuable datasets.

`seqsubmit` is an ambitious, community-driven initiative to create an **nf-core pipeline to automate submission of diverse genomic data** (samples, reads, assemblies etc.) to INSDC databases.

INSDC (International Nucleotide Sequence Database Collaboration) comprises three major repositories which mirror each other's data:

- ENA (European Nucleotide Archive)
- NCBI GenBank
- DDBJ (DNA Data Bank of Japan)

Each database has its own submission interfaces, tools, and requirements—making truly universal submission tools a significant challenge that requires expertise across multiple platforms and active community involvement.

## Goals

For this hackathon, we'll focus on creating reusable nf-core subworkflows for **metagenomic assembly and bin/MAG (Metagenome Assembled Genome) submission to ENA**. These subworkflows can be integrated with existing pipelines like [`nf-core/mag`](https://nf-co.re/mag) for automated submission.

1. Create missing nf-core modules for:

- [`assembly_uploader`](https://github.com/EBI-Metagenomics/assembly_uploader) (MGnify's tool for metagenomic assembly submission) (Difficulty: Easy/Medium)
- [`genome_uploader`](https://github.com/EBI-Metagenomics/genome_uploader) (MGnify's tool for MAG/genome submission) (Difficulty: Easy/Medium)
- [`webin-cli`](https://ena-docs.readthedocs.io/en/latest/submit/general-guide/webin-cli.html) (ENA's general submission client) (Difficulty: Easy/Medium)

2. Write a `nf-core/mag` samplesheet generator for `nf-core/seqsubmit` (Difficulty: Medium/Hard).
3. **Success**: Develop an nf-core subworkflow for uploading metagenomic assemblies with proper metadata handling (Difficulty: Hard).
4. **Absolute Success**: Develop an nf-core subworkflow for uploading MAGs and bins, taking pre-calculated quality statistics (completeness, contamination, coverage) and taxonomic assignments as input (Difficulty: Very Hard).

Let's make data submission more friendly—one subworkflow at a time! 🚀

## Resources

- [`seqsubmit` project proposal](https://github.com/nf-core/proposals/issues/79)
- [ENA Metadata Model](https://ena-docs.readthedocs.io/en/latest/submit/general-guide/metadata.html)
- [ENA Assembly Submission](https://ena-docs.readthedocs.io/en/latest/submit/assembly.html)
- [ENA Genome Submission](https://ena-docs.readthedocs.io/en/latest/submit/genome.html)
