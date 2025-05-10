---
title: Annotate all the proteins with proteinannotator!
category: pipelines
slack: "https://nfcore.slack.com/archives/C084HP11VPF"
intro_video: ""
image: "/assets/images/events/2025/hackathon-boston/proteinannotator-moar-proteinanotator-modules.jpg"
image_alt: nf-core proteinannotator hackathon
leaders:
  olgabot:
    name: Olga Botvinnik
    slack: "https://nfcore.slack.com/archives/D08JCBP586S"
---

This project will add modules to the nf-core/proteinannotator pipeline, building towards a 1.0.0 release!

## Vision

Build the best protein annotator in the world.

Protein fasta -> ??? -> Profit!

- the ??? is `nf-core/proteinannotator`
- We want to build the pipeline of choice by the people sequencing the genomes of new creatures to annotate protein fasta files with function
- Future options include using synteny of genes, but that is beyond the 1.0.0 release

## In-progress annotation tools

So far, we have these PRs in progress that we could use your help on!

- [#9](https://github.com/nf-core/proteinannotator/pull/9): [InterProScan](https://interproscan-docs.readthedocs.io/) -- Started by @olgabot, will work on during the hackathon -- **looking for contributors to update InterProScan on nf-core/modules**
- [#14](https://github.com/nf-core/proteinannotator/pull/14): Convert Fasta to Parquet files to compute amino acid composition stats using [FastaToParquet](https://github.com/heuermh/dishevelled-bio/blob/352ff5578a11a8b84755fc7b251362ee3adb847f/tools/src/main/java/org/dishevelled/bio/tools/FastaToParquet.java#L73) from [heuermh/dishevelled-bio](https://github.com/heuermh/dishevelled-bio), started by @hueuermh
- [#17](https://github.com/nf-core/proteinannotator/pull/17): UniProt's [UniFire](https://gitlab.ebi.ac.uk/uniprot-public/unifire) -- [Instructions](https://www.ebi.ac.uk/training/events/annotate-your-proteins-uniprot-functional-annotation-system-unifire/), started at the March hackathon and **looking for more contributors!**
  - This may end up needing to be written as its own _subworkflow_ because the UniFire container from EBI runs its own internal pipeline that duplicates work, e.g. it runs InterProScan internally and uses the output from that for further analysis
- [#18](https://github.com/nf-core/proteinannotator/pull/18): [DIAMOND-blastp](https://github.com/bbuchfink/diamond), started at the March hackathon and **looking for more contributors!**

## TODO annotation tools

- [proteinfold](https://nf-co.re/proteinfold/1.1.1/) + [FoldSeek](https://github.com/steineggerlab/foldseek) -- will require folding protein structures, e.g. with ESMFold2 or AlphaFold2 -- **looking for contributors!**
- ... Another tool you suggest ...?

Plus any of the below!

- BLAST/BLASTP: https://nf-co.re/modules/blast_blastp/ **(already an nf-core module!)**
- CLEAN: https://www.nature.com/articles/s41467-024-55676-y
- DHR: https://www.nature.com/articles/s41587-024-02353-6?fromPaywallRec=false
- HMMer: https://nf-co.re/modules/hmmer_hmmsearch **(already an nf-core module!)**
- MMSeqs: https://nf-co.re/modules/mmseqs_search/ **(already an nf-core module!)**
- PLMsearch: https://www.nature.com/articles/s41467-024-46808-5
- ProtTucker: https://academic.oup.com/nargab/article/4/2/lqac043/6605840
- ProtNote: https://github.com/microsoft/protnote
- REBEAN: https://www.biorxiv.org/content/10.1101/2024.12.10.627786v2.full
- S2F: https://www.nature.com/articles/s42256-021-00419-7?fromPaywallRec=false
- TM-align: https://nf-co.re/modules/mtmalign_align/ **(already an nf-core module!)**
- TM-Vec & DeepBLAST: https://www.nature.com/articles/s41587-023-01917-2?fromPaywallRec=false

_We welcome contributors of all experience levels._

## Similar pipelines

Below are pipelines that also process protein fasta files and add either functional or structural information to them, but don't have exactly the same purpose as `proteinannotator`. We will likely use their modules.

- [funcscan](https://nf-co.re/funcscan/dev/) to search (meta)genomic nucleotide data for functional protein sequences, e.g. for biosynthetic gene clusters, antimicrobial peptide genes, and antimicrobial resistance genes
- [reportho](https://nf-co.re/reportho/dev/) to compare ortholog predictions across methods
- [proteinfamilies](https://nf-co.re/proteinfamilies/dev/) to cluster protein sequences into families, and updates existing families with new sequences
- [proteinfold](https://nf-co.re/proteinfold/1.1.1/) to fold protein sequences with ESMFold, AlphaFold2
