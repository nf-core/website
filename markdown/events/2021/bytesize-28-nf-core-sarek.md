---
title: 'Bytesize 28: nf-core/sarek'
subtitle: José Fernández Navarro - BioLizard, Spain
type: talk
start_date: '2021-11-23'
start_time: '13:00 CET'
end_date: '2021-11-23'
end_time: '13:30 CET'
embed_at: 'sarek'
youtube_embed: https://youtu.be/6EIGUe5sjNo
location_url:
  - https://youtu.be/6EIGUe5sjNo
  - https://www.bilibili.com/video/BV1WR4y147PT
  - https://doi.org/10.6084/m9.figshare.17068046.v1
---

# nf-core/bytesize

Join us for a special pipeline-focussed episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 28: nf-core/sarek

This week, José Fernández Navarro ([@jfnavarro](https://github.com/jfnavarro)) will tell us all about the nf-core/sarek pipeline.

nfcore/sarek is a workflow designed to detect variants on whole genome or targeted sequencing data. Initially designed for Human, and Mouse, it can work on any species with a reference genome. Sarek can also handle tumour / normal pairs and could include additional relapses.

<details markdown="1"><summary>Video transcription</summary>
**Note: The content has been edited for reader-friendliness**

[0:01](https://youtu.be/6EIGUe5sjNo&t=1)
(host) Hi, Maxime here. I'm glad today to welcome you again for another bytesize talk. Today it will be about the pipeline Sarek and I think it's a very good opportunity that we have one of our power users who will be presenting this pipeline. In the end, it would be so boring if I were to present that. I would like to thank again the Chan Zuckerberg Initiative to help us do this bytesize talk every week. Right now, let's head over over to you, José to present us Sarek.

[0:48](https://youtu.be/6EIGUe5sjNo&t=48)
Thanks, Maxime. Hi, everyone. I would like to start by thanking the nf-core and the Sarek people for allowing me to give this talk. As Maxime said, I'm one of the power users and now also contributor. I would also like to thank the nf-core community in general. I am an active user of other pipelines and I'm always present in Slack and I'm very happy of this initiative and how the community is growing and they're setting up now the standards for pipelines. I'm really thankful. Well, this talk has been live streamed and recorded on YouTube. I like to say hi to everyone that is out there watching or will be watching later on.

[1:38](https://youtu.be/6EIGUe5sjNo&t=98)
A little bit about me, I guess this is a mandatory slide, but I'm not gonna talk too much about me because this is about Sarek, it's not about me. I have a bachelor degree in computer science, a master in computational biology and a PhD in computational biology. You want to check something about my research or my projects, you can see my GitHub or my ResearchGate page. Some facts about me, I'm born and raised in Spain. The picture there, right there is my hometown this summer. But I did most of my studies and my career in Sweden, Stockholm, and I got to meet Maxime and Phil and some nf-core people there. I'm one of the regional developers of Spatial Transcriptomics, which is a technology that you might have heard of. So I accumulated a bunch of years experience both in academia and industry. I was even group leader of a bioinformatics unit at some point, but I decided to change careers and move to industry. Now I'm senior consultant at BioLizard, which is a consultancy company that specialize in bioinformatics. But that's enough about me.

[2:51](https://youtu.be/6EIGUe5sjNo&t=171)
A first slide, that I like. Many of you probably wonder where does the name Sarek comes from. I found out that this comes from a national park in Northern Sweden. This picture is beautiful. If you have a question about why this national park and why this name, you can ask Maxime later. But Sweden is beautiful. I spent there 10 years. Sarek is a Nextflow nf-core pipeline to detect germline and somatic mutation in whole genome sequencing and whole exome sequencing data. As you know, the pipeline is written in Nextflow and it's part of the nf-core pipelines. It's heavily used, I think, and the community is growing. I see the Slack channel growing and the amount of users growing. I am one of the users, I started getting familiar with Sarek as a user, but now I'm also contributing. But the main authors are Maxime and Sylvester. I believe this project started in Sweden, Stockholm, and it was started by Maxime and Sylvester, but then a bunch of contributors have been added in the last year. I am one of them, but also Gisela, Friederike, and many others. You can see this in the GitHub page, the list of contributors. In the GitHub page the Sarek page in the nf-core. I also have to say that this work is published, or at least soon to be published. It's in peer-review, I believe, at this stage. You might also want to check the publication. It's public access.

[4:36](https://youtu.be/6EIGUe5sjNo&t=276)
Sarek, as I said before, is a Nextflow pipeline for doing variant calling in genomic data. This is a really nice illustration that is part of the Sarek documentation on the Sarek GitHub page. As you can see here, Sarek is designed for both germline and somatic workflows. This is the main workflow. This is a very high level, but Sarek does more things. But yeah, the most important thing is Sarek follows the GATK best practices, 4.0, which are the standard for pre-processing of genomic data for variant calling. Sarek, one of the cool things that I like, has a lot of tools, both for germline and somatic workflows. For example, for germline, it has HaplotypeCaller, FreeBayes, mpileup, Strelka2, Manta, TIDDIT, I don't know how to pronounce that. For somatic, obviously, Mutect2, which is very popular, FreeBayes, Strelka2, Manta, ASCAT, Control-FREEC, MSI-sensor.

[5:42](https://youtu.be/6EIGUe5sjNo&t=342)
But it does more things. This figure is a very high level illustration of the workflow, but Sarek does also quality trimming with Trim Galore, which is a wrapper around cutadapt. QC with FastQC, BamQC, a mapping step, which is essential, it's done with BWA or BWA2. As I said before, it follows GATK4 for best practices for pre-processing and marking duplicates in the realignment.

[6:09](https://youtu.be/6EIGUe5sjNo&t=369)
Also, I'd like to mention that Sarek can be used in tumor-only mode. It's a somatic mode, but you only have two more samples. Sarek has been updated, so it can work with tumor-only samples. It's also compatible with exon and targeted data, which is very nice because essentially Sarek contains almost everything. It's originally designed for mouse and human references, but technically possible to use other references. As far as I know, me personally haven't used any other reference, but there might be people out there that have used them. It would be nice if someone can say that at the end.

[6:49](https://youtu.be/6EIGUe5sjNo&t=409)
I would like to start by talking about Sarek's germline mode. When one wants to do germline variant calling, it's when we have samples that are not somatic or they might be somatic, you might want to be germline variant calling. But it's essentially to detect variants that are not in the reference genome. Let's say you have a sample and you want to see the variants that are not in the reference genome. They could be genomic variants, variants that want one inherits. In order to use Sarek in germline mode, you just need to include... because Sarek has this option to allow multiple tools to be run. If one of the tools that I provide are from the germline toolset, Sarek will do the germline workflow. These tools are HaplotypeCaller, Strelka2, FreeBayes, for SNPs and indels, Manta for structural variance, and TIDDIT for structural variance. The input for Sarek is a tab-delimited file, which is standard in the nf-core pipelines. The user needs just a tab-delimited file with information about the samples and path to the raw data. The output will be BAM files for intermediate steps and VCF files for each of the callers that are included in the run. Of course, a very nice report in what form are done with MultiQC.

[8:25](https://youtu.be/6EIGUe5sjNo&t=505)
Somatic mode, probably quite standard case for variant calling, is to detect variants using a reference genome, but also another sample, a normal sample. These will be variants that are somatic, that are acquired, that might be specific to a tumor, specific to a cell type or to an individual. For this, we need the tumor sample, but also the normal sample that can be used as a reference. How to run Sarek in somatic mode? Just as I said before, to include tools that are for the somatic workflow like Mutect2, Strelka2, FreeBayes, Manta, ASCAT, Control-FREEC for copy number variation, and MSI-sensor for MSI status. In this occasion, to run Sarek in somatic mode, the input file has to contain both the normal and the tumor sample. It has to contain this information. It has to be indicated which is which, which I will describe later, and the path to the files. The output would be BAM files, VCF files, and other files, because, for example, the copy number variation and the MSI status, they don't generate VCF files, they generate other types of files. As before, the reports in MultiQC.

[9:49](https://youtu.be/6EIGUe5sjNo&t=589)
Sarek also supports tumor-only mode. This is a case where you might want to run the workflow in somatic mode, but you don't have a normal sample, you only have a tumor. This happens sometimes. Of course, the variants are less reliable, you are more prone to have false positive, but it's something that Sarek supports. I'm actually one of the users of this mode. As before, in order to use Sarek in somatic mode, tumor-only, some of the tools that are supporting the tumor-only mode has to be included, like Mutect2, Strelka2, Manta, Control-FREEC, MSI-sensor, and a tab-delimited file as the other modes, but this time only with the tumor sample. The output will be the same depending on the tools. You always have the BAM files, intermediate files, but depending on the tools that you include in the analysis, you will get VCF files, copy number, MSIs, and the report.

[10:49](https://youtu.be/6EIGUe5sjNo&t=649)
I have mentioned before the input file, the tab-delimited file, and this is something that I see that users, especially when they start, struggle a little bit with the format of the input file. I did struggle myself. The input file format is a tab-delimited file that has a bunch of columns. The first column is the subject ID. This will be something that could be also the dataset ID or project. This will be something that will be appended to the output, to the different files, the different folders. The gender, which is used in some of the tools, and the very important column, the third column, is whether the sample is tumor or not. Zero means non-tumor, and one means tumor, and this is very important for Sarek to know whether to run the somatic mode, the germline mode or somatic mode in tumor-only. The sample ID, of course, this is important, especially, this will be appended to the output, but especially when you have multiple lanes and you have multiple FASTQ files that you want to merge for the particular sample, Sarek will take care of that.

[12:01](https://youtu.be/6EIGUe5sjNo&t=721)
For example, in this case, sample ID one has three different lanes. This will be merged before doing the processing and the path to the files, of course. As simple as this, sometimes people get confused, but this is the main required parameter for Sarek, and obviously the tools that one wants to run, as well as the reference genome.

[12:26](https://youtu.be/6EIGUe5sjNo&t=746)
One thing that I felt it was nice to explain, and for some of you this might be very obvious if you're familiar with variant calling or with Sarek, but I think it's nice to go through what Sarek does, the information that Sarek provides with different tools, because Sarek has a lot of different tools that can provide different information, so it's good to know, to have an idea. The variant calling usually involves deriving SNPs, indels, and structural variants, and so Sarek has some tools that are specifically designed for this, like Mutect2, HaplotypeCaller, Strelka2, and FreeBayes tools will compute SNPs and indels. This figure is quite representative of what this type of information is. Essentially, we have a reference genome that has a sequence, and you have a sample that you have sequenced and in the sequence of a specific region, you might see that there is a single nucleotide that is changed as compared with the reference, so this would be a SNP, but the change could also be an insertion or a deletion, so it would be indel, representing the second and the third row here, so these are variants that are detected.

[13:52](https://youtu.be/6EIGUe5sjNo&t=832)
Of course, if you're in germline mode, you only compare with the reference, but if you're in somatic mode, you would also compare with another sample, so you want to see that this change is not only happening in the reference, but it's also happening in the sample that's used as reference. Structural variants are more like bigger changes, for example, deletion or insertion, or a big portion of the DNA sequence. They could also be inversions. The tools that are providing this information in Sarek are Manta, and TIDDIT.

[14:31](https://youtu.be/6EIGUe5sjNo&t=871)
Other type of information that Sarek can provide are copy number variations or MSI status. The copy numbers are just a difference in the number of copies of specific genes, different species will have a predefined number of copies that you should have for your mother, for your mother and father, in case a human. This is a very, very common analysis in genomics to derive the copy number of a sample to see which regions have a different copy number. It could be more or less, illustrated here with red and green, and Sarek has integrated two different tools for inferring copy number status, copy number variations, region-wise, chromosome-wise, and these tools are ASCAT and Control-FREEC.

[15:21](https://youtu.be/6EIGUe5sjNo&t=921)
Another thing that is interesting, especially in cancer, is the MSI status. An MSI individual will be a hypermutated individual, they usually have defects in the DNA repair pathways, so these samples are hypermutated, they have a higher number of mutations, and this is good to know, especially in the cancer field. Sarek uses MSI-sensor, which is a tool that provides this status, which is one score per sample, and this you can see here, the percentage of MSI. Usually what one does is to put the threshold to define a sample as MSI or not.

[16:06](https://youtu.be/6EIGUe5sjNo&t=966)
Now we move to the output, the VCF format. Sure, many of you are familiar with this format. This was a format that was designed to include mutational data variants, and I personally like it, some people might not, but it's quite standard now in the field, and I think the current version is 4.3, so most of the tools I would say, the variant calling tools, they make use of this format. The format is quite flexible. It has a header that contains metadata, meta information: for example it contains the commands that were run, which is very nice to have an history of what happened in that file, how that file was generated, and also it explains all the fields, the format, the info, each explanation for every field. If you add extra things to the file, usually you have to include it in the header, so you explain what those new fields mean.

[17:07](https://youtu.be/6EIGUe5sjNo&t=1027)
Then the body. The body contains the mutational information, the first column is the chromosome, then the position in the reference where that mutation is detected. An ID, usually it's an RS number that is used to identify, to track this mutation in databases, for example. The reference of the genome, and the alter, let's say the new nucleotide in case of single nucleotide mutation. The quality of the mutation. The filter is something that many variant callers provide, they have some confidence scores, they use some times probabilistic methods to give a confidence to that variance. Some tools they provide if the filter is passed or not, which can be used to process or post-filter the vcf file.

[18:06](https://youtu.be/6EIGUe5sjNo&t=1086)
Information, this field is quite open, here the information about fields that could be added to the body, information about the body, annotation, the format also contains information about the variant, the number of reads in each allele. One thing that I like about vcf files is you can have multiple samples, so the format will detect the information of each sample that is included, and then the samples will be concatenated here. One single vcf file can contain information about multiple samples. I don't think I need to explain this here, but of course in the vcf file, different type of mutation will be represented differently, like SNPs, insertions, deletions, replacement, or structural variants, and if you want to learn more about this format, there's a really nice specification, this pdf here in the SAMtools website, and this picture was taken from this website.

[19:12](https://youtu.be/6EIGUe5sjNo&t=1152)
One of the things that Sarek does, as you could see in the workflow, is the annotation, which is something crucial in my opinion. Once you have a vcf file with the variants, you want to annotate them, which means to assign the functionality, get the extra information about the variants. Sarek currently supports two different annotation tools, VEP, which is developed by Ensembl people, I think it's variant effect predictor, stands for that, and also SnpEff, which I will describe later. Anovar is not currently supported, and I'm not sure that there are plans to include it, but I think it's a proprietary software, so that might make it difficult. But VEP and SnpEff are similar in the way the information they provide. They have differences, some information VEP can provide, some information that SnpEff can provide, and the other way around. But the most important information... I mean the annotation will create a new vcf file that contains the annotation for each of the variants. And annotations are quite important. It's the gene, the transcript, where this variant is detected, the feature type, it could even give a consequence. It could tell what consequence that variant can have upstream, the position in the genome, the position, the amino acids change, the codon. It could even predict effects. This is very useful for cancer. It could even use databases like COSMIC to tell you that variant has been detected, population information, how common is that variant in the population. This information is quite useful, especially to filter variants or to detect variants that could be potentially interesting.

[21:02](https://youtu.be/6EIGUe5sjNo&t=1262)
SnpEff is another tool, I particularly like this one, it has a really nice documentation, it's really easy to use and fast, and provides similar information: feature type, feature ID, the gene name, biotype, what kind of, is it non-coding gene, is it protein coding, RNA, the impact of the variant. This is something that is specific to SnpEff. It has different categories like high, low, moderate, which could be used to filter. They call it here putative impact, but it's the effect, the downstream effect of the body. Once again I have personally used this a lot in order to filter the data and to detect potentially interesting variants.

[21:48](https://youtu.be/6EIGUe5sjNo&t=1308)
The last step of Sarek is output. Obviously the output is not only the vcl files and the bam files, it's this really nice report. Sarek will automatically generate a web report, MultiQC which is standard in the nf-core community, this is a tool that was developed in Stockholm, Sweden by Phil, and I love it. I really like it. It's a web report, it's an HTML report that contains a lot of different stats for the different steps of the pipeline, I cannot show everything here, but yeah, essentially you get information about the mapping quality, duplicates, the type of variants that are detected, amount of variants, the region where the most variants were detected. A lot of information that is included here. It's quite useful to make an assessment of the data, an initial assessment to detect samples that might be outlayers. Samples that might be discarded. And to get a general overview of the data.

[22:53](https://youtu.be/6EIGUe5sjNo&t=1373)
One use case. I have personally used Sarek in three projects at least that I remember now. I think the biggest one was this project that I executed where I was working in Vall D`Hebron Institute of Oncology as a group leader, the bioinformatics unit. In this period we wanted to analyze 141 whole-exon sequencing samples. That's quite a lot. We used Amazon for that and Nextflow Tower. Amazon Batch for distributed computing, distributed processing, which is really nice, because you can process a lot of samples in a relatively short time. I think it took days once we had everything set up. What we did here was to process different kinds of samples, we have germline samples, we have tumor samples, some of them were solid tumors, some of them were cell-free DNA. We also have PDX samples where obviously we have to extract the human reads. But all these different samples were processed with Sarek. We used these colors that are shown here in this figure, most of the tools, and we used annotation provided by SnpEff. We then manually annotated with Anovar. We were quite happy with the result with Sarek, and this was one of the reasons why I started to contribute to it, because some of the samples here were tumor-only, and some of the features that we needed were not present in Sarek, so we just contributed. That's the best thing of open source.

[24:27](https://youtu.be/6EIGUe5sjNo&t=1467)
For this project the manuscript is under preparation, we got really nice results, and we were really happy with the short time we could process all the samples, the amount of information provided, and the quality of the variants. This is something very important when doing variant calling. It's quite tempting to just develop your own pipeline, but these tools have a lot of different settings, and it can be tricky to optimize and to run all the steps properly, and Sarek is a pipeline that is built from a community, and it has all these users. I really recommend using Sarek if you want to do genomic variant calling.

[25:11](https://youtu.be/6EIGUe5sjNo&t=1511)
Approaching the end, what is next? Sarek 3.0 is coming soon, so stay tuned. I hope that this release will happen in early 2022. I'm contributing to this release, and I'm quite happy to. One of the highlights of this new release is that it will be ported to DSL2, which is now the new syntax and the new language for Nextflow. It will include more tools. I think for now DeepVariant, I might be missing something, maybe Maxine later can elaborate here. It's been refactored and redesigned, obviously, because DSL2 forced you to do that, but now Sarek is more modular. It has a simplified design. It's easier to navigate it. It was a bit difficult before. More tests will be added. The testing will be improved, and validation tests will be included. The joint variant calling, which is something very useful, especially when using the Sarek gemline mode.

[26:14](https://youtu.be/6EIGUe5sjNo&t=1574)
Downstream analysis, this will be really nice, maybe to add the option to filter vcf files according to certain criterias, maybe merging. This is something that is now ongoing, and it might be very handy for users. Decrease the resource needed. The current version of Sarek is doing computation that are not needed in certain cases, Sarek now will be optimized, only what is needed will be used.

[26:45](https://youtu.be/6EIGUe5sjNo&t=1605)
Get involved. I recommend that. If you're a user of Sarek and you like programming and Nextflow, if you have something that is missing, I really recommend to get in contact with the community. There is a Slack channel. People are quite responsive and very helpful, very nice. There's also a GitHub page where you can see instructions on how to contribute. I really recommend that if you're using Sarek and something is missing, just at least say in the Slack channel and you might get it. You may see it happening.

[27:21](https://youtu.be/6EIGUe5sjNo&t=1641)
Last slide. The acknowledgments. As I understand, Sarek was initiated at Karolinska, at Scilifelab, which is a place that I know very well. I spent there 10 years. So thanks to these institutions for supporting. Obviously, nf-core community, I say many times, I am personally quite thankful for the effort that they are doing. Nextflow, Nextflow Tower, BioLizard, which is a company that I'm working for now and they're supporting my work now with Sarek with the 3.0 release. Vall D`Hebron Institute of Oncology, which is a place that I was working where I ran that trial that I explained before. I did some contribution to Sarek. QBiC, which I think is also a contributor, I believe Gisela. They are also supporting the development on Sarek. I think that's it.

[28:23](https://youtu.be/6EIGUe5sjNo&t=1703)
(host) Thank you very much, José, for the talk. We have a couple of questions already. I see that, okay, Laurence has a question, but you already responded to it.

(question) James has a question about, like in the current Sarek version, does Mutec2 support multi-sample calling?

(answer) I don't think it does in the current version, but definitely that's something that will be interesting. We will have a look maybe in the 3.0 or maybe like in another version. That's something we can think of.

[28:59](https://youtu.be/6EIGUe5sjNo&t=1739)
(question) Then another question from Laurence again. Do you have any idea how a SnpEff defines a putative impact?

(speaker) I cannot hear you well, sorry, can you repeat?

(question cont.) Yes, sorry. Do you know how SnpEff defines a putative impact?

(answer) I think they have some database.

(host) Yes, I think as well they have some database about that.

(answer cont.) I just had to say that you have to be careful with that because sometimes, for example, moderate variance, sometimes it can be quite interesting. How a SnpEff defines these different categories, you have to be careful because sometimes miss-sense mutations are defined as moderate. I particularly, in one project that I had, I was filtering out moderate variance and I realized that they were the most important. The one that is to use these filters, has to be careful. Essentially the filter is just including a bunch of the effects. Which effects are included in this category you can see on the website, the documentation is quite nice.

(host) Yes, I think Laurence's question was about modifier and moderate. Yes, she will double check the documentation.

[30:15](https://youtu.be/6EIGUe5sjNo&t=1815)
(host) Okay, thank you very much for the presentation. It was super clear and everything, but I might be a bit biased about that. Definitely I will steal some slides for my own presentation because that was even clearer than when I usually explain stuff. That's good.

(speaker) Yeah, thanks to you and everyone attending.

(host) Okay, then I think we are good there. Thank you very much and see you another time.

</details>
