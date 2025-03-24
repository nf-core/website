---
title: Inaugural work for nf-core/proteinannotator Pipeline
category: pipelines
slack: "https://nfcore.slack.com/archives/C084HP11VPF"
intro_video: ""
image: "https://media0.giphy.com/media/v1.Y2lkPTc5MGI3NjExa3JhcDNlcWR4N3czMTI4ajhvcXNvNXdmZWt3enlnc216M3lveXZqOCZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/xUOxfjsW9fWPqEWouI/giphy.gif"
image_alt: nf-core proteinannotator hack
leaders:
  olgabot:
    name: Olga Botvinnik
    slack: "https://nfcore.slack.com/archives/D08JCBP586S"
---

This project aims to begin work on the nf-core/proteinannotator pipeline.

## Vision

Build the best protein annotator in the world.

Protein fasta -> ??? -> Profit!

- the ??? is `nf-core/proteinannotator`
- We want to build the pipeline of choice by the people sequencing the genomes of new creatures to annotate protein fasta files with function
- Future options include using synteny of genes, but that is beyond the 1.0.0 release

BEFORE WRITING ANY CODE, we will first draw out the metromap for

## Similar pipelines

Below are pipelines that also process protein fasta files and add either functional or structural information to them, but don't have exactly the same purpose as `proteinannotator`. We will likely use their modules.

- [funcscan](https://nf-co.re/funcscan/dev/) to search (meta)genomic nucleotide data for functional protein sequences, e.g. for biosynthetic gene clusters, antimicrobial peptide genes, and antimicrobial resistance genes
- [reportho](https://nf-co.re/reportho/dev/) to compare ortholog predictions across methods
- [proteinfamilies](https://nf-co.re/proteinfamilies/dev/) to cluster protein sequences into families, and updates existing families with new sequences
- [proteinfold](https://nf-co.re/proteinfold/1.1.1/) to fold protein sequences with ESMFold, AlphaFold2

## Annotation Tools to Include

Please contribute more tools! This is just a starting point.

- [DIAMOND-blastp](https://github.com/bbuchfink/diamond)
- [InterProScan](https://interproscan-docs.readthedocs.io/)
- UniProt's [UniFire](https://gitlab.ebi.ac.uk/uniprot-public/unifire) -- [Instructions](https://www.ebi.ac.uk/training/events/annotate-your-proteins-uniprot-functional-annotation-system-unifire/)
- [FoldSeek](https://github.com/steineggerlab/foldseek) -- will require folding protein structures, e.g. with ESMFold2 or AlphaFold2

_We welcome contributors of all experience levels._
