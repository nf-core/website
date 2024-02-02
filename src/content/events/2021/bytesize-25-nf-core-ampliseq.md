---
title: 'Bytesize 25: nf-core/ampliseq'
subtitle: Daniel Straub - QBic, University of Tuebingen, Germany
type: talk
startDate: '2021-10-26'
startTime: '13:00+02:00'
endDate: '2021-10-26'
endTime: '13:30+02:00'
embedAt: 'ampliseq'
youtubeEmbed: https://youtu.be/a0VOEeAvETs
locationURL:
  - https://youtu.be/a0VOEeAvETs
  - https://www.bilibili.com/video/BV1B44y1e7MM
  - https://doi.org/10.6084/m9.figshare.16871008.v1
---

This week, Daniel Straub ([@d4straub](https://github.com/d4straub/)) will tell us all about the nf-core/ampliseq pipeline.

nfcore/ampliseq is a bioinformatics analysis pipeline used for amplicon sequencing, supporting denoising of any amplicon and, currently, taxonomic assignment of 16S, ITS and 18S amplicons. Supported is paired-end Illumina or single-end Illumina, PacBio and IonTorrent data. Default is the analysis of 16S rRNA gene amplicons sequenced paired-end with Illumina.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://youtu.be/a0VOEeAvETs&t=1)
(host) Hi, everyone. I hope you're all ramped up for the nf-core hackathon starting tomorrow, but it's time for bytesize now, and I'd like to begin by thanking the Chan Zuckerberg Initiative for funding nf-core events. I'd also like to thank those of you from the community that have actively been helping shape our outreach events recently and for making our material accessible. We're back this week with a bytesize talk that will be focused on pipelines, and we have Daniel Straub, who is based at QBiC at the University of Tübingen in Germany, talking to us about the nf-core/ampliseq pipeline. If you have any questions for Daniel, please either unmute yourself at the end of the talk or use the chat function, and I will read out your questions. Thanks very much for taking the time to present at bytesize, Daniel. Over to you.

[0:46](https://youtu.be/a0VOEeAvETs&t=46)
Thanks, Renuka. As she already mentioned, I'm talking about nf-core/ampliseq today. To give an introduction, microbial communities are very important in our everyday life. For example, we as humans have more bacterial cells than human cells. In agriculture, bacteria or microbial communities are also very important. They can be either beneficial or pathogenic, so it's also quite important here to know what's going on in order to improve. Also in the environment, microbial communities play major roles, like they are driving the global elemental cycles, meaning that the elements are turned over and recycled and can be used again by nature to build whatever it can build. For example, here on the right side is a bacteria that is in the iron cycle, global cycle, or here on the lower figure, we also have a microbial community which is involved in removing arsenic from groundwater in Asia, where we also investigated the microbial community.

[2:06](https://youtu.be/a0VOEeAvETs&t=126)
There are different ways to investigate microbial communities or mixed communities, mixed organisms at nf-core. We have, for example, a pipeline called nf-core/bacass, which is for genome assembly of isolated bacteria. From a microbial community, you can isolate a bacteria, extract their DNA, and then assemble its genome, or from the whole community, can total community RNA be extracted, and then the meta-transcriptome can be investigated with nf-core/rnaseq. The DNA can be also extracted from the total microbial community to assemble the metagenome with the nf-core/mag pipeline, or a fragment of the community DNA can be amplified by a PCR, polymerase chain reaction, and produce amplicons that are typically done for ribosomal RNA genes, rRNA genes, and that kind of data can be analyzed with nf-core/ampliseq that I'm talking today about.

[3:23](https://youtu.be/a0VOEeAvETs&t=203)
Not only ribosomal RNA gene fragments can be amplified, but any other kind of gene or fragment that you might be interested in. However, ribosomal RNA gene amplification is typically used to investigate microbial communities and is the default of nf-core/ampliseq, so I will focus now on that. Ribosomes are essential for the life we know it. Ribosomes are translating the RNA in cells to peptides, to proteins. Ribosomes contain proteins but also non-translated RNA, and that is the ribosomal RNA, which is important for its structure. Prokaryotes have a specific set of ribosomal RNA, including the 16S rRNA, which is used mostly for investigating the taxonomic composition of microbial communities. Eukaryotes have different kinds of ribosomal RNA, and that is also used to investigate communities of those.

[4:29](https://youtu.be/a0VOEeAvETs&t=269)
The 16S rRNA gene that encodes for that non-translated 16S rRNA, that is highly conserved, but it has regions that are variable and other regions that are very similar. These variable regions of the 16S rRNA gene allow the discrimination of many bacterial taxa. Having conserved regions means that we can use PCR primers for the PCR to amplify that variable region. The problem of that is that there is an error accumulation during library preparation and sequencing. You can see that sample sequences here as four different dots in different sizes corresponding to their abundance have errors, meaning different kind of sequences when they are sequenced. The sequencing reads have errors. The traditional way of correcting this is producing so-called operational taxonomic units, OTUs, and those OTUs are lumping together all sequences that are similar to a specific threshold level. For example, here it would then produce only two OTUs out of that four different original sample sequences. There are other methods which make so-called amplicon sequencing variants, ASVs, and those methods try to correct the errors of the amplicons back to the original sequences that are found in the samples.

[6:13](https://youtu.be/a0VOEeAvETs&t=373)
First thing I did is to compare those kind of methods. What you can see here on the y-axis is the relative abundance of sequences in a reference data set. The reference data set consists of 27 bacterial strains, which has 35 expected amplicon sequences because bacterial genomes can have several copies of that gene. Here you can see the reference, meaning this is how it should look like. There are high abundant sequences, but also very low abundant sequences. The color of the dots show if the sequence is a perfect match in green to the reference, in blue one off, so very similar, or in black cross a sequence that is not expected. You can see here that the OTU-based methods, here four columns for four different methods or settings, have quite a lot of those unexpected sequences. They also, in this region here, do not find all the ones that we would expect, the green dots here. The ASV-based methods, here two different methods, do show a much reduced number of false positives. However, some of those sequences are also not detected for this method here, that is very strict.

[7:50](https://youtu.be/a0VOEeAvETs&t=470)
I did that not only with one sample, but with three different samples and with here the OTU-based methods and the ASV-based methods, not only for precision and sensitivity of finding the original sample sequences again, but also for the taxonomic classification and for a diversity index, the Shannon diversity. Here the green color shows the best methods and the red ones, the worst ones. You can see that the ASV-based methods perform much better in the amplicon sequencing analyses than the OTU-based methods.

[8:31](https://youtu.be/a0VOEeAvETs&t=511)
Coming from that, I wanted to produce a pipeline that is doing exactly what I need and to use that kind of knowledge that I have accumulated by the benchmarking. I wanted to produce a pipeline to analyze all kinds of amplicon data with a wide range of input types, reference taxonomies, and downstream analyses. We came up in 2018 then with the nf-core/ampliseq pipeline, in the end of 2018. We released version 1.0.0 and did some feature and maintenance updates. We also then published a paper about it and continued afterwards in 2021 quite a lot with adding more features to the pipeline. Now it is at a pretty good stage how I imagined it in the beginning.

[9:31](https://youtu.be/a0VOEeAvETs&t=581)
What is happening in nf-core/ampliseq? And so this is the whole picture. This is a bit complicated, so we will go step by step through. First of all, to start the pipeline, you would need such a command here. In red are the ones that are specific to that pipeline. First there are different kinds of inputs that can be used. Either direct input of a folder containing all your reads or a sample sheet in a TSV. That is what I describe here. The TSV has to contain at least two columns with sample ID and forward reads and two optional columns with reverse reads and sequencing run, because different sequencing runs of a sequencing machine has to be treated separately by the pipeline. That sample sheet is then added with the input parameter. Using that command will then end up in primer trimming, QC, quality filtering, and then also the inference of this amplicon sequencing variants. For this, of course, you also need to add your primer sequences because it is very important what kind of sequence was amplified.

[10:57](https://youtu.be/a0VOEeAvETs&t=657)
Depending on what kind of sequencing data from what kind of sequencing technology you have, you have to also add additional parameters. For Illumina paired-end sequencing, just default, you don't need to add anything else. For single-end, the `--single_end` parameter, for PacBio and for IonTorrent technologies, you have to specify that. In special cases, Illumina paired-end ITS, so fungi community analyses, you also should already set a parameter so that all the downstream analyses will have sensible settings.

[11:38](https://youtu.be/a0VOEeAvETs&t=698)
What is the output of those very basic analyses? In the results folder, in our subfolder "dada[2]" will be that file, and that file here contains already most of the information that you will need. It contains for each of the different samples, the sequences that it found and also the quantification, so it simply counts how many of the sequences in there will originate from that original 16S sequence in your sample. There is also a handy overall summary for all the read numbers that were for each of the samples processed and then ended up in the table and if those numbers do not look as expected, then this is a very good starting point to troubleshoot.

[12:34](https://youtu.be/a0VOEeAvETs&t=754)
The next step of the pipeline is to classify the sequences that we have now produced, meaning that we give it a name, that we give it a taxonomic name. This is by default done with the DADA2 software, but we also have an alternative way to do this with QIIME2 which is also a very popular program in the area. By default it is using the DADA algorithm with the SILVA138 taxonomy. So we have here a range of reference taxonomy databases that can be used and that is the SILVA, the RDP, the GTDB database that are all for bacteria, for 16S rRNA gene amplicon analyses; the UNITE database for fungi analyzes ITS, internal transcribed spacers, and the PR2 database for eukaryotes with the 18S rRNA gene amplicons.

[13:41](https://youtu.be/a0VOEeAvETs&t=821)
You can also taxonomically classify any kind of sequences that you have produced by using `--input` and provide your FASTA, which needs to have then a .fasta extension. Following, there will be taxonomic filtering and by default only mitochondria and chloroplasts will be filtered out from the following tables, because those are typically off-targets of the 16S rRNA gene amplicon sequencing, that should be removed. The output of this taxonomy classification is again a table, a TSV, that will list all the sequences that were discovered and also add here the taxonomic levels and their classifications with confidence value.

[14:37](https://youtu.be/a0VOEeAvETs&t=877)
Finally there is also downstream analysis, but this downstream analysis is only provided when the optional parameter metadata is pointing to a tab-separated metadata sheet. Then automatically appropriate metadata columns will be chosen and used for visualization of a bar plot, as you can see here. Those are HTML files that are then interactive so you can choose colors, you can choose different kinds of taxonomic levels. I chose here a very high one just to show you, for example here with the test data set, that we are typically running here with triplicates, that you can see already patterns of different kind of data.

[15:27](https://youtu.be/a0VOEeAvETs&t=927)
Additionally there will be a differential abundance analysis with a program called ANCOM and that produces for example such a volcano plot, also interactive. This red dot here which is identified as a significant different abundance between the different treatments would be a _Burkholderiales_ bacteria.

[15:52](https://youtu.be/a0VOEeAvETs&t=952)
Then also alpha and beta diversity indices are produced. Alpha diversity is a measure of each of the samples, how diverse the sequences in there are. You can see that for this test data set, it's a bit small, the samples originating from groundwater, have the lowest Shannon diversity index but the sediment samples here have the highest one. To all of these comparisons will be of course also statistical analyses provided in the output of the pipeline. Beta diversity plots show the difference in PCOA plots. What you can see here is that the replicates for the different kinds of sample sources (groundwater, river water, sediment and soil in our example data set) are nicely clustering and well separated from the others.

[16:55](https://youtu.be/a0VOEeAvETs&t=1015)
Additionally there will be quality control figures like alpha rarefaction curves, which show you if you have a sequence with sufficient output for your samples. As long as those curves will flatten out you know that you have sufficient sequencing data produced for that sample. All those results can be also found on the nf-core website in [`/ampliseq/results`](https://nf-co.re/ampliseq/results) and you can have a look what is produced, how does it look like and there is also extensive documentation what each of those plots mean and what you can use it for. Those plots, of course, cannot show all kinds of complicated experimental setups, but they are very helpful for getting a first view on the data and are also a sort of quality control that what you are producing here makes sense. With this I want to finish and I want to thank you for your attention, my colleagues at QBiC , who also were involved in the development of the pipeline, nf-core, especially Daniel Lundin, helped a lot and the University of Tübingen Microbiology and Geomicrobiology who produced a lot of data, that I was then able to analyze and with this pipeline and produce it as good as possible.

[18:33](https://youtu.be/a0VOEeAvETs&t=1113)
(host) Thanks very much for that very clear and insightful talk Daniel. Are there any questions? If anyone in the audience has questions you can unmute yourself or even add questions to the chat and I can read them out. Okay I don't see anything coming up so that was extremely clear so I guess that's why people don't have questions but if anything crops up anywhere, you know where to reach Daniel. You can message him and continue the discussion on the bytesize channel and thanks again for joining. We will be back next week with another pipeline focused talk.

[19:23](https://youtu.be/a0VOEeAvETs&t=1163)
(question) We have one comment here okay there's a question. This person thanks you for the presentation and asks what type of license is associated with this kind of tool, for example the ampliseq pipeline.

(answer) It is a CCPY license [comment: ampliseq is actually under an MIT licence], like I think all the nf-core pipelines have, there are no restrictions whatsoever.

(host) I hope that answers your question. Yes thank you.

[20:02](https://youtu.be/a0VOEeAvETs&t=1202)
(host) As I said we will be back with another pipeline focused bytesize talk next week where we will have Payam Emami, who will be presenting the metaboigniter pipeline. But before that I hope to be able to see some of you in sunny Gathertown, starting tomorrow at our hackathon and take care everyone thanks again for joining.

</details>
