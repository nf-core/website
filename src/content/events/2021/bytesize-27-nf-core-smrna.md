---
title: 'Bytesize 27: nf-core/smrnaseq'
subtitle: Lorena Pantano - NextRNA Therapeutics, USA
type: talk
start_date: '2021-11-09'
start_time: '13:00+01:00'
end_date: '2021-11-09'
end_time: '13:30+01:00'
embed_at: 'smrnaseq'
youtube_embed: https://youtu.be/4YLQ2VwpCJE
location_url:
  - https://youtu.be/4YLQ2VwpCJE
  - https://www.bilibili.com/video/BV17U4y1M7bQ
  - https://doi.org/10.6084/m9.figshare.16964392.v1
---

# nf-core/bytesize

Join us for a special pipeline-focussed episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 25: nf-core/smrnaseq

This week, Lorena Pantano ([@lpantano](https://github.com/lpantano/)) will tell us all about the nf-core/smrnaseq pipeline.

nf-core/smrnaseq is a bioinformatics best-practice analysis pipeline used for small RNA sequencing data.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://youtu.be/4YLQ2VwpCJE&t=1)
(host) Hi, thanks for joining us for another bytesize talk that's focused on pipelines. I'd like to begin by thanking the Chan Zuckerberg initiative for funding nf-core events and our outreach initiatives. I'd also like to extend a special thanks to the volunteers that helped us organize the recent hackathon two weeks ago. Now, we've recently received funding from the Chan Zuckerberg Initiative that is focused on boosting diversity and inclusion. Evan mentioned a few full-time paid positions for developer advocates for this purpose in the announcements channel on Slack. Please reach out to us for more information and keep an eye out for the official advertisement of the positions and spread the word when they're out. Now, without further ado, I'd like to introduce you to our speaker today. We're joined by Lorena Pantano from NextRNA Therapeutics in the US, who will be presenting the nf-core smRNA-Seq, so small RNA-Seq pipeline. If you have any questions for our speaker today, please unmute yourself at the end of the talk or use the chat function and I will read out your questions. Lorena, we appreciate how early in the morning it is for you. Thank you very, very much for making the time to present for us at bytesize. Over to you.

[1:14](https://youtu.be/4YLQ2VwpCJE&t=74)
Thank you. Thank you for having me in mind to share my contribution to the community. I'm very happy to be here and it's a sunny day. It's pretty awesome, the view. I'm going to talk about this pipeline, the small RNA-Seq pipeline. What I'm going to do is to give a scientific background, before starting to talk more about the analysis. I always think about the small RNA as small non-coding RNA. There are different kinds of small non-coding RNA that some of them are handled by the pipeline, others currently not, or maybe there are better pipelines for that.

[2:00](https://youtu.be/4YLQ2VwpCJE&t=120)
When we talk about small, it's good to give a size range. It could go even to a little bit shorter than 20, maybe 18 to 400. They are involved in different biological functions. We have gene regulation where we have microRNAs that are very well known. Together with microRNAs, we have endogenous siRNA and piRNAs, that are a little bit longer. We have others that are more in the range of 100, 200 that are involved in post-transcriptional modification like polyA or RNA nucleotide modification and others with similar size in the translational machinery like transfer RNA, or ribosomic RNA. Then there are others that we aren't clear of the function but as well are notated and are part of this group of small molecules. I'm going to go a little bit over the ones that the pipeline can handle, although not all of them right now have tools that will quantify them properly, but it's something that we want to work towards.

[3:18](https://youtu.be/4YLQ2VwpCJE&t=198)
In general when we talk about these very small RNA molecules we are talking about sizes between 18 to 33. Normally the biogenesis can come from double strand RNA or hairpin structures and it will be processed probably by Dicer, but not always, an Ago machinery, to perform some gene silencing. But it's not always like that but it's a good summary of what it could be. The microRNAs, these are the most well studied small RNA, the precursor normally is a big primary precursor that contains several hairpin structures and this long precursor is cleaved just around this hairpin structure to form the pre-microRNA. This still is in the nucleus and then it will go to the cytoplasm to get the final mature microRNA thanks to Dicer. At the same time of this process it will find what RNA it will target by nucleotide complementarity. There is a lot of research done in that this silencing can go through can go to mRNA degradation or translational repression.

[4:50](https://youtu.be/4YLQ2VwpCJE&t=290)
The one thing that maybe not a lot of people know is about the variability or variations that you find in microRNAs. It's very well known that microRNAs are not just one single RNA molecule, for everything we know there are variations, and it's not only single nucleotide variations. There are as well variation at the beginning and at the end of the sequence. These are called isomiRs and they were first described, I think, in 2008 and there is a lot of research on that. This modification is important, it will change the function. They will correlate to different stages in the cells, to diseases, so it's important to understand that when you want to analyze this you want to differentiate at least the type of variations that you have.

[5:50](https://youtu.be/4YLQ2VwpCJE&t=350)
Endogenous siRNAs have been described a lot in C. elegans and Drosophila. They are involved in regulation and they have a similar biogenesis than microRNAs. They can come from hairpin or double-strand RNA even in cis or trans, but at the end they will go to the cytoplasm and will perform a similar regulation. From here is where the field of using siRNA exogenous sRNA to knock down genes comes.

[6:26](https://youtu.be/4YLQ2VwpCJE&t=386)
Finally but not least, although there are a few other categories, piRNAs. They appear later and are a little bit longer from 26 to 30 nucleotide long. They have some modifications at the end of the sequence and were discovered the first time by regulating the transposons in the genome. When everything gets de-methylated, many of the transposons are active and this is a way to regulate that. They are created in this feedback loop. A small RNA that will have this size and it will amplify, it will then bind to the DNA to induce methylation or direct the degradation of the RNA from the transposons. It got a lot of interest in the past years.

[7:28](https://youtu.be/4YLQ2VwpCJE&t=448)
Other types of small RNA are coming from tRNAs, actually. There is a lot of investigation in here. They have been used a lot as biomarkers in cancer or different diseases. The function is still unclear, although there is a lot more information. But the idea of the biogenesis, when you do small RNA-seq and analyze the data, if you quantify these RNAs, you will find that there are fragments, small RNAs that are belonging to one side or another, or the top or the bottom. There is high abundancy, so there is no way that this is degradation because it's a conserved part of the whole molecule and they have been correlated to many diseases already.

[8:28](https://youtu.be/4YLQ2VwpCJE&t=508)
Then in the same way that there are fragments coming from tRNA there are fragments coming from other nuclear RNAs. In the same way this is still less clear. There is some conservation of the processing of these RNAs so this is where people keep quantifying and see what could be the function. As you see, many of them you already know, but this is not only about microRNAs, there is a lot of types of RNA and ideally, in a comprehensive pipeline you want to try to quantify and to characterize all of them. It's very hard to put all that together so I will go through the analysis, how it will look like and what is integrated in the nf-core pipeline and what is more to come.

[9:24](https://youtu.be/4YLQ2VwpCJE&t=564)
There is a lot of people working on that at different times and I hope that more people will come, because as you see it's not that simple a case to analyze this data. One thing that is recurrent is, because small is not a super quantitative word, that you will think that, oh, you have 200 nucleotides sequencing data and you think it's small and you will try to use this pipeline. I wouldn't recommend that. I think that if your library is enriched in reads for mainly more than 50, the best pipeline to use is RNA-Seq. Maybe there should be another pipeline for this, but all the tools that this small RNA-Seq pipeline has integrated is based on that you are going to find a lot of duplication of the reads, because your molecules are shorter than the read length you normally have.

[10:30](https://youtu.be/4YLQ2VwpCJE&t=630)
I will say that this is a good yes-or-no decision you're making. If you have a mixed one. I have done a lot that I have analyzed the data with this pipeline and then the untrimmed reads I will put in the RNA-Seq because sometimes you have a very heterogeneous library. You can do that, it's fine.

[10:54](https://youtu.be/4YLQ2VwpCJE&t=654)
The main challenges of the pipeline is isomiR detection because there are a lot of changes in a very small region. Right now we have some strategy, but ideally we will integrate more. There are smaller RNAs that come from places that are in multiple parts of the genome, so it's difficult to handle that and depending on the library you may have difficulties to differentiate what is degradation or functional molecules. Also, when you go to non-model organisms things get a little bit difficult. Ideally a perfect pipeline should try to address all this.

[11:35](https://youtu.be/4YLQ2VwpCJE&t=695)
The first step will be trimming, which is related to the library prep. There are different library preps. Normally there is variation of the two main ones, which are TruSeq, where you attach adapters to your RNA or cDNA molecule in this case, or you have a NextFlex library, where you have some degenerated bases at the end and at the beginning of the cDNA so you need to take that in account. The pipeline is compatible with a few protocols but probably not all the combinations. That is another issue we keep having all the time. We try to adapt to to all the possible combinations.

[12:23](https://youtu.be/4YLQ2VwpCJE&t=743)
The pipeline is divided in processing and QC, detection, annotation of non-microRNAs and trying to discover the novel ones. In processing and QC we use TrimGalore. We have FastQC integrated in there. MirTrace is a... I will give an example of this, but this is a good QC tool and of course MultiQC to put everything together. I would like to integrate qualimap to give you an idea of any other small RNA that you have, although MirTrace helped with that already. For detection and annotation we have bowtie and samtools as an easy way to know how much is mapping to mature, precursor or genome. Mirdeep2 will give you a quantification of microRNAs and mirtop will focus on isomiRs. I really want to integrate MINTmapper or any other tool that handles tRNA, there are some options. Because as I say, they have become very important. For the novel detection we have mirdeep2. It will do the job , but there are other tools like seqcluster or protac for piRNAs that will try to handle other kinds of RNAs.

[13:38](https://youtu.be/4YLQ2VwpCJE&t=818)
An example of MirTrace. It has a lot of different QC metrics but this one I like because it gives you an idea of the diversity of your sample and it will quantify what belongs to microRNA, to ribosomal RNA, to tRNA, to PCR artifacts or primer dimers. It gives you quickly an idea of what samples are better or if you have anything that you need to remove. For the quantification, as I say, after trimming and MirTrace you will go to mirdeep2 and bowtie and samtools and what you are going to have is a matrix of microRNA and isomiRs from these tools. This is what you will use for differential expression. Mirtop is compatible with the isomiR R package in bioconductor and it gives you a lot of tools to understand isomiRs, like a map of the different isomiR abundances in different samples and sometimes you see differences. It is helping you to decide where to focus.

[14:55](https://youtu.be/4YLQ2VwpCJE&t=895)
Finally this is the part that needs more work right now. Mirdeep2 is running different samples. A better study will be to collapse all the samples into one, to run mirdeep2, so you can really detect everything at once and then separate the quantification. Right now it is running on the different samples. Here is where seqcluster and other tools, that are focused on characterizing other non-coding RNAs, will be good to be added at some point. From this you will get the matrix and you can work with visualization and differential expression as you wish. An example of mirdeep2. It will give you the secondary structure of the RNA where the mature is. This visualization is pretty cool and it is good to go over that to see if you have any new microRNA in your samples.

[16:02](https://youtu.be/4YLQ2VwpCJE&t=962)
Just as a summary, the current pipeline has a very comprehensive QC with MirTrace. The microRNA isomiRs are very well quantified. But it needs a couple of more tools to be a better pipeline. The novel microRNA detection and quantification with mirdeep2 is there, so that is good. We are migrating to the second version of Nextflow. It has been massive work for many people and we are almost there. We are working with some final details, so hopefully we can get it done in a few weeks. Mainly it is the versioning and adding some more QC metrics. My idea is, because I have this data in my current job, that it gives me an excuse to work on this. I would really like to add other tools to make it a better pipeline and if anybody wants to help they are very welcome. For this or any other pipeline I am happy to help in any way as well.

[17:23](https://youtu.be/4YLQ2VwpCJE&t=1043)
As I say this is a work of many. I couldn't put everyone here but I really appreciate everybody. I really appreciate the Nextflow and nf-core community and I think that is helping a lot of people: beginners, experts, everybody. So I really like the vision and the core values. I want to thank contributors as well, people enabled in this.

</details>
