---
title: 'Bytesize: nf-core/airrflow'
subtitle: Gisela Gabernet, QBiC, University of Tübingen
type: talk
startDate: '2022-10-25'
startTime: '13:00+02:00'
endDate: '2022-10-25'
endTime: '13:30+02:00'
youtube_embed: https://www.youtube.com/watch?v=CrJgxVRVlqY
locationURL:
  - https://doi.org/10.6084/m9.figshare.21407004.v1
  - https://www.youtube.com/watch?v=CrJgxVRVlqY
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: nf-core/airrflow

This week, Gisela Gabernet ([@ggabernet](https://github.com/ggabernet)) will talk about the newest developments in the nf-core/airrflow pipeline.
nf-core/airrflow is a bioinformatics best-practice pipeline to analyze B-cell or T-cell bulk repertoire sequencing data. It makes use of the [Immcantation](https://immcantation.readthedocs.io/) toolset and requires as input targeted amplicon sequencing data of the V, D, J and C regions of the B/T-cell receptor with multiplex PCR or 5' RACE protocol.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=1)
Hello everyone, my name is Franziska Bonath and I'm very happy that Gisela is with us today from the University of Tübingen and she is giving us an overview of what nf-core/airrflow can do.

[0:13](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=13)
Thanks Franziska for the nice introduction.
I'll present nf-core/airrflow. First of all I will start with defining what's the airr in airrflow. Airr stands for adaptive immune receptor repertoire and that's the collection of membrane proteins that are found on the surface of B-cells, in which case they are called BCR, and on the surface of T-cells, in which case they are called TCRs or T-cell receptors. BCRs or B-cell receptors in their secreted form are also called antibodies, which is a term that we are all more familiar with.

[0:52](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=52)
The main function of these receptors is to be able to recognize foreign antigens that are inside the human body, which can come for example from pathogens such as viruses or bacteria, and to elicit an immune response against them. To be able to recognize so many different antigens from different pathogens, these receptors have to have a variety of different sequences. It's estimated that in the human body at any time point, there's 10 billion to 100 billion different receptor sequences of BCRs and TCRs.

[1:27](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=87)
Airr sequencing is about getting the individual sequences of these B-cell and T-cell receptors and that has a variety of applications, which can vary from determining the immune state of one individual at a specific point, studying immune related diseases, guiding vaccine development or guiding cancer immunotherapy. I'm going to talk about a bit more detail on how this diversity of the TCRs and BCRs is generated.

[1:58](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=118)
In the human genome there are different V-gene segments, D-gene segments, J-gene segments and C-gene segments, that can be combined to form these TCR and BCR receptors. In case of humans there's four V-gene segments, 23 D-gene segments and six J-gene segments. What happens in the B-cells and T-cells is that each one segment of each kind is combined to form a productive TCR or BCR receptor. That happens at the DNA level in a process called somatic recombination. It alters the genome of the cells and generates this productive TCR and BCR receptors.
The combination of the different gene segments doesn't happen like lego blocks, just attaching them next to each other, but rather in a cut-and-paste procedure. The cutting position is not always exactly at the same spot and there can also be some nucleotides incorporated into these junction regions, so that there's extra variability that comes from this step, that is not genome encoded.

[3:12](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=192)
In the case of the B-cells or in B-cell receptor, there's even another process that generates more diversity of these receptors, which happens upon antigen stimulation. That's what happened to all of us, for example, when we were first in contact with the coronavirus or the coronavirus vaccine. Some B-cells in the body were able to recognize the antigens in the coronavirus and they were stimulated and underwent clonal expansion that is generating a lot of children cells that belong to the same B-cell clone. This happens in a manner that not all of these children cells have the exact same BCR receptor sequences. The process called somatic hypermutation introduces mutations in these VDJ segments, so that each of the children cells has a slightly different receptor sequence. This allows to generate B-cell receptors and therefore antibodies, that have even higher affinity to the original antigen.

[4:17](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=257)
You've seen that there's many different processes that contribute to the diversity of this TCR and BCR sequences, including somatic recombination, variable junction length and in the case of the BCRs also somatic hypermutation. This means that theoretically there could be more than 10 to the power of 14 possible BCR sequences.

[4:42](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=282)
How are they sequenced? Most protocols use what's called amplicon sequencing, which is the targeted amplification of this gene locus and that can be done via different protocols including multiplex PCR, where one is providing primers for all the kinds of different sequences that can be there. Another technique is five prime race amplification protocols, which are pretty common. This sequencing protocols can incorporate what's called unique molecular identifiers which allow correcting for sequencing errors down the line and errors introduced by the PCR amplification process.
Typically for sequencing, MiSeq sequencers are used, because they allow for longer read lengths, that covers the complete VDJ and beginning of the C region.

[5:39](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=339)
How is the bioinformatic analysis done for these kind of sequences? This typically does not happen like a traditional RNA-Seq analysis and that's because mapping to a reference genome is challenging in this case, due to the high diversity of these BCR and TCR sequences. Also, it's a highly repetitive genome region with all of these VD and J gene segments there. For this kind of analysis, specific tools are used and the reads are aligned to specific reference data, also for these BCR and TCR receptors.

[6:15](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=375)
I always say luckily for us, when wanting to write a pipeline to do this kind of analysis, there's already plenty of tools out there that can analyze this data. One of the better well-known ones is the Immcantation framework that is developed by the Kleinstein lab in Yale. It's an open-source toolset to analyze AIRRseq data from beginning to end. There's a whole community of users already using this framework and here you can also get the details, in case you want to have more information.
Thanks to developing this pipeline as part of the nf-core community, we gained visibility and we quickly found quite some collaborators to develop the pipeline. I want to mention that this is really a community-based development effort! Susanna Marquez from the Immcantation lab joined early in the beginning and also David Ladd from Monash University helped by adding some features. Whereas initially, Alex Peltzer,
when he was still at QBiC with us, Simon Heumos and myself were developing the airrflow pipeline.

[7:26](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=446)
Now to the details of the pipeline, those are the main pipeline steps. The pipeline can process both bulk AIRR sequencing data and single-cell sequencing data. When starting from bulk, there's first a step of quality control of the sequencing reads and sequence assembly and afterwards there's a process, where the reads are aligned to the references with IgBLAST. The reference data is typically employed from the IMGT consortia, which provides reference data for BCR and TCR. At this step already assembled data can also be provided and for single-cell data typically we start at this step.
Afterwards, there's a step for clonal analysis, which identifies which of the sequences of the BCR belong to which B-cell clone. It assigns the individual BCR sequences and TCR sequences to their specific clone and in the case of the B-cells it can also perform lineage reconstruction of the whole B-cell clone. Finally there's a step for reporting, doing repertoire analysis and reporting, including QC reports by MultiQC. That's the general steps of the pipeline, but now if we look a bit more detailed there is a ton more processes - individual processes that are part of the pipeline. I'm going to explain them a bit in more detail now.

[9:01](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=541)
Starting for the QC and sequence assembly, the pipeline supports different sequencing protocols, including multiplex PCR, in which the users have to provide the V and C-primer sequences that were used for the amplification, or in the case of 5' RACE, providing the C-primer and the linker sequences for amplification. Both protocols are supported with and without UMI barcodes. The barcodes can also be provided in different configurations. Starting from the raw sequencing data, a sample sheet needs to be provided that contains sample information and the individual FASTQ files for all samples. Depending on if the sequencing protocol includes UMI barcodes or not, there will be some processes but they all start with quality control of the reads with FastQC, filtering the sequences by quality threshold, masking the primer sequences and if a UMI-based protocol is used, then a consensus is built from all the sequences that have the same UMI barcode. This way, it also allows to correct for the errors as I mentioned before.
There is an extra procedure, that is employed whenever it's estimated that the length of the UMI barcodes will not be sufficient to cover all of the diversity of the sample and that is bypassed by first clustering all the sequences by similarity and annotating the cluster ID and then two different sequences with the same UMI barcode can also be distinguished this way.

[10:51](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=651)
After building the consensus, all the sequences that contain the same UMI barcode are collapsed and the count, the number of sequences with the same UMI barcode, is annotated. Other metadata is also annotated as well as the count of duplicate sequences with different UMI barcodes that were collapsed. This can be useful for filtering. After this step there comes the VDJ assignment and filtering step and here it's also possible
to start with already assembled sequences, that can be provided with a sample sheet and FASTA files. Typically single cell sequencing data processing starts at this step, because the pipeline supports directly the output from the tool 10x Genomics CellRanger multi, which provides, while incorporating also TCR and BCR sequencing in the 10x Genomics sequencing procedure, the output of that tool, which is the AIRR rearrangement
table, that contains all of the sequences there. This can also be directly provided to the pipeline at that step. What that step does, it aligns the sequences to the IMGT reference and then assigns what exact V, D and J segments were used there. In the case of the single cell data there is an option. This gene reassignment step is optional.

[12:29](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=749)
After alignment to the IMGT reference there is a number of quality filtering steps that are performed. First it's checked that the locus matches exactly the V_calls. The V(D)J segment assignment, that there are a minimum of 200 informative positions, maximum 10% and nucleotides, that the sequences that are determined are productive, that the junction region is a multiple of three amino acids. There's also a possibility of removing chimeric reads, detecting contamination across samples and finally collapsing duplicate sequences if there are any. Plenty of quality filters as part of the pipeline.

[13:18](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=978)
The next step is clonal analysis. In that case hierarchical clustering is used based on the hamming distances between sequences. The pipeline is also able to auto-detect a hamming distance threshold that can be used to determine, which sequences are part of the same clone, or are part of different clones. There is also a step for lineage reconstruction of the clonal lineage trees. Recently the pipeline offers to use the EnchantR tool, which is developed by Susanna Marquez in Immcantation. That tool provides calls to other Immcantation tools and also nice reports for each of these steps. I invite you to check them out here.

[14:05](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=845)
Finally, there's a repertoire analysis step. An R Markdown report is provided, that summarises the repertoire analysis results for all samples. Of notice is that a custom R Markdown report can also be provided, in case that the user wants to change some things in this report. It's also possible to provide an R Markdown file. Other reports of this reporting analysis steps is the MultiQC QC report for all samples from the grid quality control reports.

[14:43](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=883)
Here we will see an example of this repertoire analysis report. First there is a summary of all the samples used for the analysis. Clonal abundance and clonal diversity are reported, together with vision usage. Finally all of the tools that are used as part of the pipeline and their citations are noted here, to make it really easy for users of the pipeline to also cite the original tools, that are being used.

[15:12](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=912)
As you know, all documentation for nf-core/airrflow can be found on the nf-core website, so check it out if you want to use it. There are also some example results of the pipeline when the full tests are run on AWS.

[15:30](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=930)
What's next for this pipeline? Stay tuned for a new release that comes really soon - we hope this week or the next. It will includes more quality control and reporting as part of the EnchantR tool, as I have mentioned, by Susanna Marquez, as well as and code refactoring using subworkflows.

[15:50](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=950)
At this point I would like to thank all the contributors to the pipeline: Simon Heumos and myself at QBiC, Alex Peltzer who was initially at QBiC but now is at Böhringer Ingelheim, Susanna Marquez at the Kleinstein Lab in Yale, David Ladd at Monash University and some collaborations also at the University of Tübingen, Christoph Ruschil and Markus Kowaric.
If you have any questions don't hesitate to join us, join the airrflow channel on nf-core Slack and if you have any questions related to the Immcantation tools, you also have the contact emails to contact them directly.

[16:26](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=986)
(host) Thank you very much. Everyone can now unmute themselves if there are any questions. Maybe I start with one, I'm curious. In what format does airrflow expect UMIs to be provided? Does it have to be in a separate file or should it be in read one or in read two?
(answer) It supports all kinds of these configurations that you have mentioned. You can provide them. It depends on your library design, where the UMI barcodes are located. Sometimes they are part of the R1 and R2 reads and sometimes they are part of the index reads, it can be provided in any way.
There are some parameters in the pipeline where you can specify where the UMI barcode is located: R1 reads, R2 reads or index files, so everything is supported in that case.

[17:23](https://www.youtube.com/watch?v=CrJgxVRVlqY&t=1043)
(host) Thank you.
Are there any questions from the audience? It doesn't seem so. Then I would like to thank of course Gisela but also the Chan Zuckerberg Initiative for funding the bytesize talks, and as usual, if you have any questions go to Slack at nf-core/airrflow and ask questions there. Thanks again.

</details>
