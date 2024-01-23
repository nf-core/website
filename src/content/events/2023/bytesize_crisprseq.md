---
title: 'Bytesize: nf-core/crisprseq'
subtitle: Júlia Mir (QBiC, Tübingen University) and Marta Sanvincente (Pompeu Fabra University, Barcelona)
type: talk
startDate: '2023-02-14'
startTime: '13:00+01:00'
endDate: '2023-02-14'
endTime: '13:30+01:00'
embedAt: 'crisprseq'
youtube_embed: https://www.youtube.com/watch?v=x_eFQW0nNvo
locationURL:
  - https://doi.org/10.6084/m9.figshare.22100045.v1
  - https://www.youtube.com/watch?v=x_eFQW0nNvo
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: nf-core/crisprseq

This week Júlia Mir ([@mirpedrol](https://github.com/mirpedrol)) and Marta Sanvincente ([@msanvicente](https://github.com/msanvicente)) present the newly released nf-core pipeline nf-core/crisprseq.
Nf-core/crisprseq is a bioinformatics best-practice analysis pipeline for the analysis of CRISPR edited next generation sequencing (NGS) data. It allows the evaluation of the quality of gene editing experiments using targeted NGS data.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=x_eFQW0nNvo&t=1)
Hello everyone to this week's bytesize talk. I'm very happy to welcome today Julia from QBiC in Tübingen and Marta from UPF in Barcelona. They're going to talk about another new pipeline that was released just a week ago called crisprseq. Off to you.

[0:22](https://www.youtube.com/watch?v=x_eFQW0nNvo&t=22)
Thank you. Thanks for the introduction. We'll present nf-core/crisprseq, which is a pipeline for the analysis of CRISPR experiments. I would like to start by an introduction to what CRISPR is, because I'm sure you've heard that word before, but maybe you don't remember exactly what it is. CRISPR comes from bacteria and the system is repurposed to do gene editing. It consists of a protein that we call Cas, and this protein can cut DNA, creating double strand breaks. It's coupled to a single guide RNA, which is a short sequence of RNA, which is complementary to the DNA region that you want to cut. This way we can have directed cuts.

[1:18](https://www.youtube.com/watch?v=x_eFQW0nNvo&t=78)
When we have this double strand break in a cell, there are usually two mechanisms of repair. The most common one is this one that we call non-homologous end joining. That's the cell that goes and tries to repair this double strand break, and this can produce some insertions or deletions, which can result in the disruption of the gene, and then this can cause a gene knockout. Then there's a different way, which is called homology-directed repair, which consists on having a template that we can provide and the repair is made based on that template. Like this we can introduce new fragments of DNA and possible gene knock-ins. Apart from these two mechanisms, there's also this microhomology-mediated end joining, which is very similar to the non-homologous, but it happens when there are two small regions of homology surrounding the cut, and these can recombine, so we can get a bigger deletion. More recently, there are these other two technologies called base editing and prime editing, which are done not by a double strand break, but only with a nick. Those are more precise because they can produce base substitutions of only one base.

[3:04](https://www.youtube.com/watch?v=x_eFQW0nNvo&t=184)
That's the overview of all these CRISPR-Cas experiments that we can have. Apart from that, we can also have CRISPR screens, which consist of a library of different gRNAs targeting lots of different genes, and then we can perform a screening. Finally, if we couple with a CAS protein that's inactive and doesn't cut the DNA, it only affects the expression of the gene, we call this CRISPR activation or CRISPR interference. Our pipeline, crisprseq, can analyze gene knockouts, knock-ins, and also base editing or prime editing experiments. This pipeline is based on a pipeline called CRISPR-Analyzer, which Marta developed, so she'll explain more about it.

[4:04](https://www.youtube.com/watch?v=x_eFQW0nNvo&t=244)
As Julia has already said, this first release of nf-core/crisprseq pipeline is based on CRISPR-Analytics. Currently, we just have the core of CRISPR-Analytics in crisprseq, which I will show you here. These are the core steps of that pipeline. The first steps are quality pre-processing of the sequencing reads, where different steps are done to remove low-quality reads, and also in the case that we have paired-end sequencing reads, the reads are merged. Then the alignment against the amplicon reference is done, and after that, there is a process where each indel and substitution that could be caused by these genome editing tools are quantified. Finally, some plots and tables are done to allow us to visualize the results. In the next slide, what I want to show you is other optional steps that CRISPR analytics have that are not currently in crisprseq, but we hope that we will be able to add it in the following versions.

[5:33](https://www.youtube.com/watch?v=x_eFQW0nNvo&t=333)
Just briefly, the first optional step that we have is the ability of using unimolecular identifiers to cluster the sequences, and through that clustering processes we can remove sequencing and amplification biases, as well as correct sequencing errors. We also have implemented a step that allows us to identify the amplicon reference, looking for it in a genome of reference. Then in the bottom part, you have two other steps that has allowed us to increase the precision of our pipeline. The first one is the size bias correction, in which we have implemented a simple model where we used spike-in controls of different sizes and known abundance that were used to model biases related to the amplification, with the sequence size, since longer deletions will lead to shorter sequences that will be amplified more times than longer ones. Then if we also sequence mock samples or a negative control, we can use this sample to subtract errors that can be also represented in our treated samples.

[7:17](https://www.youtube.com/watch?v=x_eFQW0nNvo&t=437)
You can choose the alignment that you want to use in the alignment step, but we have been exploring with simulated data sets the performance of different alignments together with the following part of quantifying the different edits. What we have done is to optimize the parameters of minimap to achieve better results related to the identification of the indels produced by the double strand break repair mechanism. In the following slide, we have just some examples of CRISPR-Analytics being used to analyze a bunch of samples. We have analyzed samples from three different cell lines that were edited with CRISPR-Cas9. In the first plot we see that the main pattern observed among all the insertions that have been found are homology insertions, which means that the same insertion that is in the cleavage site has been also added in this repair process. This happens with higher frequency when the nucleotide that we have free is an adenine or a thymine. As in the other two plots, what we have been exploring is the precise outcomes, which are those outcomes that are shown in a higher frequency. In that case, we also observe that among these precise outcomes, we have these homology insertions, and also we have some deletions of a cytosine when this cleavage site is surrounded by cytosines, and also we can see some micro-homology patterns that have lead to longer deletions that have also a higher representation in these samples.

[9:43](https://www.youtube.com/watch?v=x_eFQW0nNvo&t=583)
CRISPR-Analytics has been benchmarked using several datasets. We have used real data as well as simulated data, and we also created a ground truth dataset to be able to also have this dataset for the benchmarking. This ground truth dataset was generated by several collaborators, which had different subsets of reads, and they were classifying the indels that were found in the reads as indels produced by errors or indels produced by genome editing tools. Finally, these subsets have been used to calculate the percentage of addition of those samples, and we have extrapolated this percentage to calculate the distance between the percentages reported by different tools and the real distance or the established percentage of addition with this ground truth dataset. From this, we just want to highlight that our tool has good precision without relying on the addition windows. Most of the tools use a window where the edited indels have to take place to avoid reporting false positive events.

[11:38](https://www.youtube.com/watch?v=x_eFQW0nNvo&t=698)
How you can use nf-core/crisprseq? Basically, you can use the typical nextflow command, where you provide an input sample sheet, the output directory of the profile that you want to run the pipeline with, and then we also have this one single parameter to provide the aligner, by default we're using minimap but you can also choose between bwa or bowtie2, and the reason why we don't have more parameters is because most of them are provided with the sample sheet because they are dependent on the sample.

[12:23](https://www.youtube.com/watch?v=x_eFQW0nNvo&t=743)
That's how our sample sheet looks. You have the sample name, fastq_1 and fastq_2. If you have only single-end sequencing data you can only provide fastq_1. Then you provide the reference sequence, here it has been shortened for space issues. This reference is the reference the reads will be aligned to, so it's the region where you directed your cut. You also provide the proto-spacer which is the guide RNA that you used in your experiments to direct the cut. Finally, in case that you performed a homology directed repair experiment, you can also provide the template.

[13:11](https://www.youtube.com/watch?v=x_eFQW0nNvo&t=791)
That's the structure of the output folder, I won't go to all the directories in detail, but you will find all the outputs of all the tools used for pre-processing like to join paired-end reads. Also the quality filtering steps because we remove sequencing adapters, we remove low quality reads and mask low quality bases also, and then you also have the output of the alignment and finally the most important folder, which is this one called `cigar`. It's called like that because we parse the edits using the cigar field from the mapping. In these directories you will find some tables and summary tables of the edits and also plots.

[14:15](https://www.youtube.com/watch?v=x_eFQW0nNvo&t=855)
This is an example of the output plots. We report data quality, meaning that you will have a percentage of reads that have good quality, also the ones that were aligned against the reference. We also report the number of reads that were wild type or the ones that contained indels and from these indels we also classify by filter of quality and if they are located in the expected pick on the cut site and if they are above the sequencing error rate or not. Finally there's also classification between insertions and deletions, if there are insertions produced by our template and also if these indels are in frame or out of frame because the ones that will be out of frame are the more probable to disrupt a gene function and produce a gene knockout.

[15:33](https://www.youtube.com/watch?v=x_eFQW0nNvo&t=933)
Finally these further steps as Marta already commented that they are already implemented in CRISPR-Analytics and we will add them to crisprseq. This unimolecular clustering step to reduce PCR duplicates or sequencing biases, because usually in the sequencing methodology, shorter reads are sequenced more often but this doesn't mean that you have this particular long deletion more represented in your sample so we can correct with UMIs. Then also the automatic identification of a reference and some noise handling and finally also thinking already about version 2 of crisprseq about the idea that we will be able to analyze other kinds of CRISPR experiments such as CRISPR screening. If you have any doubt or want to work with us, Laurence is currently implementing this part of the analysis so you can join the Slack channel and there ask and that's it. Feel free to join this channel, test out the pipeline and see if there's something that you would like to also include. Also check the repository. Thank you.

[17:08](https://www.youtube.com/watch?v=x_eFQW0nNvo&t=1028)
(host) Thank you very much, that was a very nice talk. Are there any questions in the audience to either Julia or Marta? You can unmute yourself now if you want to, or you can write a question in a chat and I will read it out. There currently seems to be no questions but may I ask one?

(question) So one of the biggest issues that I know of CRISPR is off-target effects, but as I understand you're mapping to fairly short references, just a target. Is there any way how we could figure out if there are off-target effects with this pipeline or is there anything planned in the future?

(answer) This pipeline it's not really thought to be able to detect off-target effects. The experimental steps are based on amplification of your expected target and then you sequence with Illumina or other next-generation sequencing platforms. What you can do is, for example, if you use some prediction of which are the the targets that are more susceptible to be off-targets you can also amplify these off-targets and make the same analysis and see if there are indels in that regions.

(question cont.) But you would need to know what to look for then obviously.

(answer cont.) Yeah, we would have to add guide-seq or other analysis pipelines that you use for another experimental protocol to do the computational analysis. It's something that can be implemented in further steps.

(host) Thank you.

[19:11](https://www.youtube.com/watch?v=x_eFQW0nNvo&t=1151)
(host) Are there any other questions in the audience? If not I would like to thank you two for this great talk. I also would like to thank the Chan Zuckerberg Initiative who is funding our bytesize talks. If anyone has more questions to both of you, you can always go to slack and check either in the channel for crisprseq or you can also ask in the bytesize channel and I'm pretty sure the two will have a look at your question. Thank you very much.

(speaker) Thank you.

</details>
