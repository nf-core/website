---
title: 'Bytesize 31: nf-core/dualrnaseq'
subtitle: Regan Hayward - Helmholtz Institute for RNA-based Infection Research , Germany
type: talk
start_date: '2022-02-01'
start_time: '13:00+01:00'
end_date: '2022-02-01'
end_time: '13:30+01:00'
embed_at: 'dualrnaseq'
youtube_embed: https://www.youtube.com/watch?v=-J3Cbetk8Pk
location_url:
  - https://www.youtube.com/watch?v=-J3Cbetk8Pk
  - https://doi.org/10.6084/m9.figshare.19927178.v1
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 31: nf-core/dualrnaseq

This week, Regan Hayward ([@reganhayward](https://github.com/reganhayward/)) will tell us all about the nf-core/dualrnaseq pipeline.

nf-core/dualrnaseq is specifically used for the analysis of Dual RNA-seq data, interrogating host-pathogen interactions through simultaneous RNA-seq.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=-J3Cbetk8Pk&t=1)
(host) Hi, everyone. Thanks for joining us. And as usual, I'd like to begin by thanking our funders, the Chan Zuckerberg Initiative for supporting outreach events by nf-core. Just a minor detail before I start today's session. This talk will be recorded, and it is being recorded at the moment. The video will be uploaded to our YouTube playlist. I will be sharing the link on our website and on Slack. So don't worry, if you've missed it, you can catch up later. I'm delighted to tell you today that we're joined by Regan Hayward, who is based at the Helmholtz Center for RNA-based Infection Research in Germany. He will be presenting the nf-core dualrnaseq pipeline. It's a pipeline that is used to interrogate host pathogen interactions through simultaneous RNA-seq. There will be time for questions at the end of Regan's talk. And you can either use the chat function at any time today or unmute yourselves at the end of the talk and ask them directly. Thanks for joining us today, Regan. I'd like to hand over to you now.

[0:56](https://www.youtube.com/watch?v=-J3Cbetk8Pk&t=56)
Thanks for the introduction. Right, so the structure of my talk is going to be pretty similar to a lot of the other bytesize talks. I will be talking about a little bit of background, first of all, some mention of the pipeline, and then some future directions as well. I think it's important to start with what is dualrnaseq. It's from the RNA sequencing, it's simultaneously capturing, in this instance, a bacterial pathogen infecting a host cell. And through bioinformatic means, we're able to assign the bacterial reads to the bacterial transcriptome and the host reads to the host transcriptome.

[1:43](https://www.youtube.com/watch?v=-J3Cbetk8Pk&t=103)
There are some challenges with read assignment with dual RNA-seq data sets. I'm going to run through a few of these on the left-hand side. For example, here we have a read, which is being assigned to gene A. We're quite confident that we can happily assign that read gene A here. When the read overhangs the gene slightly, we're still pretty confident that we can say that this read belongs to gene A. When we have multiple annotations overlapping and the read occurs within this overlap, this is a little more challenging. Perhaps we want to just say that this read is being assigned to gene A, or is it ambiguous or just a little bit more complicated? Or is it ambiguous or do we want to assign a proportion of the read to gene A and another to gene B? When we read multi-maps, multiple genes, what do we want to do in this instance? Do we want to count it at all? Do we want to count a proportion or do we want to say it's gene A and gene B?

[2:50](https://www.youtube.com/watch?v=-J3Cbetk8Pk&t=170)
These multi-mapping reads are a bit of a challenge, especially when we can concatenate the genomes for dual RNA-Seq studies. We've got an intra-species challenges. This illustration I show on the left-hand side is both for host and pathogen reads. Depending on the infection ratio, generally we have a much lower proportion of bacterial reads in the sample. It becomes really important to try and assign as many bacterial reads as possible and as accurately as possible.

[3:20](https://www.youtube.com/watch?v=-J3Cbetk8Pk&t=200)
I'll spend some time in the next couple of slides talking about the bacterial transcriptome architecture, which I think a lot of people probably aren't aware of. A lot of people are probably more aware of the host side, eukaryotic side, and how splicing occurs and turn to splicing events. With the bacterial transcriptome architecture, bacterial genes are grouped into operons. On the right, we have a monocystronic operon with a single gene and a polycystronic operon with multiple genes inside. The operons are generally flanked by 5- and 3-prime untranslated regions. And within an operon, genes are co-transcribed into an mRNA transcript.

[4:12](https://www.youtube.com/watch?v=-J3Cbetk8Pk&t=252)
In Illumina-based sequencing, the transcripts are fragmented into reads, and the reads are assigned to genetic features. This example here, the red reads aren't being assigned, but the blue-assigned colored reads are being assigned to a particular gene. This brings about a few challenges. For instance, many of the bacterial annotations aren't actually complete. If you look outside the model organisms, such as E. coli, salmonella, and bacillus, a lot of the annotations don't include any of the UTR regions or small RNAs and even complete genomes. A lot of the bacterial species would be in just contigs or scaffolds.

[5:02](https://www.youtube.com/watch?v=-J3Cbetk8Pk&t=302)
In addition, a lot of the bacterial genomes contain a number of highly repetitive bacterial sequences as well, that can be difficult to assign reads to. I'll go into that a little further. We've worked out a uniqueness score per gene in each bacteria. It's a Kmer-based approach, and we're looking at each gene. We assign a number of Kmers to each gene, and depending on the uniqueness of those Kmers it can be seen in other genes, and then we can assign a uniqueness score per gene, which is each of these dots. If the Kmer is a unique, that means the gene gets a uniqueness score of one. If the Kmers appear in another gene, that means that gene will become a duplicate. You get varying levels of uniqueness per gene, which is indicated by the color, red being a duplicate, and gray not. We've set a cut off at about 50 percent saying that anything below is considered to be repetitive.

[6:14](https://www.youtube.com/watch?v=-J3Cbetk8Pk&t=374)
For chlamydia, which is a gram negative obligate intracellular bacteria, most of the genes are quite unique. If we expand this to other bacteria such as Mycobacterium leprie, Streptococcus pneumoniae, Salmonella, Tifurium, and Orientia, to Tutsugamshi, we can see quite a difference. If we look at the contrasting Orientia, which is a gram negative obligate intracellular, causing scub typhus, if you're interested. It contains a lot of repetitive elements as you can visually see here. We'll list them to a table form, where the duplicates are the ones in red, and the repetitive is anything below a cut off of 50 percent uniqueness. Chlamydia has a total of eight genes that would be challenging to assign reads to, and Orientia for example has over 1600, I think which is about 60 or 62 percent of the genome. Something to keep in mind.

[7:21](https://www.youtube.com/watch?v=-J3Cbetk8Pk&t=441)
Host genomes also contain a lot of repetitive elements as well, which I think a lot of people are probably quite aware of. For example, the mouse host genome contains about 45 percent repetitive elements, the human genome between 50 and 70 percent, and in some of your RNA-Seq studies, if you're looking into, we're looking at 70 percent. In some of your RNA-Seq studies, if you're looking into different plants, for example, the maize genome has over 80 percent of transposable elements.

[7:54](https://www.youtube.com/watch?v=-J3Cbetk8Pk&t=474)
Most of your RNA-Seq experiments typically will use genome-based approaches, such as STAR and maybe featurecounts or HT-Seq, something like that. In this pipeline we're introducing a transcript-based approach within the pipeline. We're using Salmon for this. Salmon has two modes. It's got an alignment-based mode, that uses an existing tool, let's say STAR, to align the reads, and then it'll use that bamfold to quantify. The second mode is selective alignment, which does a pseudo-alignment steepening quantification, so it's contained in the Salmon itself. And one of the advantages of using Salmon is it uses this expectation maximization algorithm. It's not just salmon, I should say, other software such as Callisto, ExpressRCM, they use this algorithm as well. It's going to assist in assigning some of these multi-mapping reads and also reads to some of these repetitive sequences, and does this through an iterative process.

[9:08](https://www.youtube.com/watch?v=-J3Cbetk8Pk&t=548)
I'll give you a really basic example. The first step would be assigning all of the uniquely mapped reads, so ones that have a really high confidence, and then it would go through an iterative process assigning the remaining reads. For example, after the first step, maybe gene A has zero reads and gene B has 100 reads, and there's a read that could be assigned to gene A or gene B, and would have a much higher probability of this read being assigned to gene B. And so using Salmon, we've found through a lot of benchmarking, is adventageous for assigning reads from dual RNAseq data.

[9:47](https://www.youtube.com/watch?v=-J3Cbetk8Pk&t=587)
We'll talk about the pipeline a little bit now. As input we have Illumina-based sequencing reads, and host and pathogen genome, and annotation. We've got FastQC, and therefore some quality control steps, and also adapter removal and read trimming through BBDuck and CutAdapt, and some pre-processing steps of merging the host and pathogen genomes and references to create this chimeric reference that we use. For parallel read mapping quantification steps, we've got a more traditional genome-based approach, which is using STAR and HDSeq in this instance, and then we have our two transcriptome-based approaches using STAR and Salmon confinement.

[10:48](https://www.youtube.com/watch?v=-J3Cbetk8Pk&t=648)
For the traditional genome-based approach, you can just apply a host genome and annotation and pathogen genome annotation, but for the transcriptome-based approaches, you will need to supply a host transcriptome, and the pathogen transcriptome is created automatically in the pipeline. The pipeline output separate host and pathogen features in various reports and plots to include correlation plots for both the host and pathogen samples per condition. Also a proportion, you also get number of reads showing the number of uniquely mapped host reads, uniquely mapped pathogen reads, multi-mapped host and pathogen reads, cross-mapped reads, so cross-mapped between species, un-mapped reads, and trimmed reads. You also get the biotype breakdown per sample as well, and depending on which method you use, or if you use a combination, you get this output for each method.

[11:55](https://www.youtube.com/watch?v=-J3Cbetk8Pk&t=715)
Status of the pipeline. There's a number of performance improvements I need to include. I need to update to the latest template as well, which we'll be doing shortly. And I'd like to include some additional outputs, it's graphical outputs and some data output as well. An example of that is WIG files that are separated by the host and pathogen. I'd like to migrate to DSL2. And probably one of the questions I get asked the most when I talk about this pipeline is, is there support for all dualrnaseq datasets? And to answer to that at the moment is no, it's just bacterial, because based on this transcriptome architecture I spoke about earlier. We're considering support for viral host pathogen datasets. If anyone's interested in this particular bit, I'd be curious to talk to them about some of the features that they would like to see. I can try and include those in my next update. I'd like to thank everyone, HIRI, the nf-core team and community, and the Salmon development team as well, they help. Thank you.

</details>
