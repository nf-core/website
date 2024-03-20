---
title: 'Bytesize: nf-core/chipseq'
subtitle: Jose Espinosa-Carrasco - Center for Genomic Regulation, Barcelona
type: talk
startDate: '2022-07-26'
startTime: '13:00+02:00'
endDate: '2022-07-26'
endTime: '13:30+02:00'
embedAt: 'chipseq'
youtubeEmbed: https://www.youtube.com/watch?v=59I-4wB1z4c
locations:
  - name: Online
    links:
      - https://www.youtube.com/watch?v=59I-4wB1z4c
      - https://doi.org/10.6084/m9.figshare.20418171.v1
---

This week, Jose Espinosa-Carrasco ([@JoseEspinosa](https://github.com/JoseEspinosa)) will talk about new developments in the nf-core/chipseq pipeline.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=59I-4wB1z4c&t=1)
I'm jose Espinosa-Carrasco from the comparative bioinformatics group, the CRG, and today I will introduce you to the nf-core chipseq pipeline.

[0:14](https://www.youtube.com/watch?v=59I-4wB1z4c&t=14)
Some background about me. I'm currently a postdoctoral research fellow in Notredame's lab at the CRG, the comparative bioinformatics group. As some of you may know this is a group where nextflow was created by Paolo di Tomaso. It's my second time I am there, so when nextflow was developed I was already there. I put these two things here because we are actively contributing to BovReg, which is a consortium to annotate the genome of the cow. My boss likes to say it's encode for cows. This is under the broader umbrella of the EuroFAANG which is a functional annotation of animal genomes, which aim is to annotate animal genomes. We are essentially using nf-core pipelines for this. Both in BovReg and in EuroFAANG. I'm also a core member of nf-core.

[1:39](https://www.youtube.com/watch?v=59I-4wB1z4c&t=99)
A little bit of background about chipseq. Probably all of you know about this, but what we want to obtain when we do chipseq experiments is this figs which show us where our transcription factors are binding in the genome, or the instance modifications. This is normally how the experimental procedure is done: There is a crosslink between the transcription factors (that are proteins) and the DNA in the place that they are sitting. This normally done with formaldehyde. After this there is a sonication procedure to get rid of the rest of the dna, then there is an immunoprecipitation step, which is what gives the name to the technique. This way we take the transcription factors that we are interested in and then the dna is purificated and the library is prepared. I'm not a wetlab guy, so probably many of you can explain this better than myself.

[2:56](https://www.youtube.com/watch?v=59I-4wB1z4c&t=176)
As i said before, this is the kind of thing that we obtain after we have run our nf-core chipseq pipeline, or other pipelines. Some figures from the nf-core chipseq pipeline. I was looking at it yesterday, and in terms of stars it is the third most popular pipeline. It has not been updated for a long time, as we will discuss today. It's a quite popular pipeline and also a quite used pipeline. It was originally developed by Chuan Wang and Phil Ewels and then it was modified to be in nf-core by Harshil Patel

[3:58](https://www.youtube.com/watch?v=59I-4wB1z4c&t=238)
This timeline shows i think if i'm not mistaken... I was looking at this yesterday, because i thought that chipseq was one of the first pipelines to be released in nf-core. It's not the first but it's among the 10 first ones. So as you can see here, it was first released in june 2019 and this is the release cycle of the chipseq pipeline itself. It was first released, as i said, during 2019 and then it was updated in November 2019, and version 1.2 was in July to 2020. These are two minor releases, this means that since this point it has not been any real big update on the pipeline.

[4:59](https://www.youtube.com/watch?v=59I-4wB1z4c&t=298)
We are working on the development of the DSL2 version of the pipeline. Actually most of the things that I will discuss today can be applied both to the DSL2 and the DSL1 pipeline, but if they can only be applied by to one of the versions, it will be to the DSL2 version, even if it's not yet the stable version. We have been dying to release the pipeline for a long time so we are approaching Sarek or even worse than Sarek and we have not released the version 2.0 (although we are very very close to it).

[5:41](https://www.youtube.com/watch?v=59I-4wB1z4c&t=341)
Here is the pipeline overview. It starts with your .fastq files, an input and a spreadsheet, that i will discuss during the presentation. There are some quality control processes like fastqc here also the adapters are removed with Trimgalore and then the alignments are performed. In the version 1.2 the only aligner that was available is bwa and now in the new version these other aligners will become available with Bowtie2, STAR and Chromap. After the alignment some aligner statistics are calculated using Samtools. Then there are these other processes, that are shown here. The replicates, if there are replicates, are merged using Picard, then duplicates are marked also using Picard. There is some quality control at the alignment level using preseq and Picard and also the bam files are filtered. The duplicates that we have marked in the previous step and also the blacklisted regions are removed.There are some regions in the genome that are difficult to align. These are known and and these regions are removed with Samtools. After all these steps some other alignment so this after all these procedures some alignment statistics are calculated. Then here we have some of the DL analysis that the pipeline performs, so we produce this fingerprint plot and the redistribution profiles.
This could also be seen as quality control approach because you can see the distribution of the profiles of the peaks, for instance, binding to the DNA. Also this strand cross correlation peaks procedure is is run with phantompeakqualtools and .bigwig files are produced with the peaks, so that they can be used downsteam and for the visualization. We also call broad or narrow peaks. The pipeline allows to have these two modes. Normally narrow peaks are called for transcription factors and broad peaks are are called for instance modifications, because the regions tend to be much wider. Then we run homer, to see where the annotation peaks that are produced are found relative to the genomic features, for instance genes. There is also this process with MACS2, which is to call consensus peaks across a given IP. We run Subread/featureCounts to have the number of reads that we found by peak.
We run DESeq only for quality control, so only the pca is used. In the previous versions, differential expression analysis was done but we agreed that this downstream processes should not be in the in the main pipeline. That's why they they have been removed.

[10:19](https://www.youtube.com/watch?v=59I-4wB1z4c&t=619)
Here i'm listing the main dsl2 updates in future. Of course the pipeline has been ported to the DSL2 syntax. This means that all models that were not yet available in modules have been implemented. We also need to implement some new modules for tools, we need several in one process, but you probably are familiar with this. More specific to the pipeline, the files containing the blacklisted regions, that i mentioned before, have been updated. Qe have included these new aligners: BWA is the default one, but now you can choose from this list. Actually I'm not entirely sure if Chromap is working as expected, probably here I will need the help of someone more familiar with this aligner. The effective genome size logic has been refactored. This is a parameter that s needed for MACS2, to annotate peaks and we have changed the logic. The input sample sheet format has been modified and, as i mentioned before, the differential expression analysis has been removed from the consensus peak comparisons of the pipeline. Of course we have fixed some bugs.

[11:45](https://www.youtube.com/watch?v=59I-4wB1z4c&t=705)
This is to show you something about these blacklist regions. As you see here the issue is closed, because this has been already implemented. But I want to throw a warning, that if you are still using version 1.2.2, you probably need to update the blacklist in the case that you are using one of the genomes where these lists are available, using the blacklist parameter. If you are using the development version you don't need to care about this.

[12:24](https://www.youtube.com/watch?v=59I-4wB1z4c&t=744)
As i mentioned before, MACS2 needs this effective genome size. This is the macs_gsize parameter that is encoded in the pipeline and we have included it now in the iGenomes configuration, the --macs_gsize size for the corresponding read length. We have calculated this based on this link here. If the genome is in the iGenome file and you provide the read length, it will be automatically taken from these maps. That's why we need this new --read_length parameter. If this is not the case, in the same way that we calculated these values, the pipeline will calculate the values for your genome, using the cage unique merge model, that has been implemented.

[13:29](https://www.youtube.com/watch?v=59I-4wB1z4c&t=809)
This is how the input looks like. You have the sample, fastq1, fastq2, antibody, control. We have seen this several times in similar formats during the bytesize talks. As you can see here, we have the sample. These samples will be merged, so for instance these two samples will be merged. Everything that is before this rep 1 and rep 2 is identical. This will tell the pipeline to merge the samples. If you have a single N rep, like in this case, you just provide the file here. If you have a paired-end, you will have to provide the second file here. This is the IP and this is the control. The control, as you can see here, is then listed here and of course has not this control field.
To run the pipeline, once you have this sample sheet, you just need these parameters. This is to run the test_full. It is taken from this link, which is what the test_full is using to run the test_full data test. You provide the genome and now you have to provide the read length so that we can take the value of the macs_gsize parameter from the map in the iGenomes. With this command you will be able to run the pipeline. In this case i used the development branch, because, as i mentioned, it's almost in production and it should be quite safe to use with the new features it has. Probably it will be released before the end of summer.

[15:23](https://www.youtube.com/watch?v=59I-4wB1z4c&t=923)
There are more parameters that you can take a look at in the parameters docs, to parameterize your run of the pipeline. Please take a look there and if you have any questions, just drop us a line in slack.

[15:38](https://www.youtube.com/watch?v=59I-4wB1z4c&t=938)
This is something that pops up many times in slack and that's why I put it here. You need controls for running the chipseq pipeline. We know that there are experiments that are old and maybe they did not have controls, but the pipeline currently is designed to be used with controls. There is a hack (that's why i put this this answer from Harshil here): you can use the ATACseq pipeline if you don't have controls, using these parameterization. In principle it should work, but ideally you should use controls. If if you are designing your experiments, you should have your controls.

[16:28](https://www.youtube.com/watch?v=59I-4wB1z4c&t=988)
This is the outputs of the pipeline with the command line that i previously mentioned. This is available also on the website. You can go there and see all the results with the full test data set. This corresponds still to the version 1.2.2, but hopefully soon it will be updated.

[16:55](https://www.youtube.com/watch?v=59I-4wB1z4c&t=1015)
We already have plans for future releases, because we wanted to get the 2.0 out and also we wanted it to be quite similar to the version 1.2.2, so that we can identify any bug or any problem that we have. Trom there we can start growing the version 2.0, if there are features that are needed by the community. These are the things that are planned for version 2.1. It will be to include the metro map. As you have seen i have done this schematical before which was not very nice. But I didn't have time to have a look at james' talk and to create it. We also would like to add the irreproducible discovery rate that is used to check consistency between replicates. It's kind of a standard because it's the measure that was used by encode. Of course we are open to ideas and if you find a bug, please tell us.

[18:02](https://www.youtube.com/watch?v=59I-4wB1z4c&t=1082)
With this I'm done. We have now a a summer break in terms of bite size talks, until the 13th of September, but probably Franziska knows better than me. If you have any questions tell me, and that's all.

(Fran) Thank you very much. I have now enabled for everyone to unmute themselves. If there are any questions you can do so and just ask them directly or put them in the chat.

[18:36](https://www.youtube.com/watch?v=59I-4wB1z4c&t=1116)
(Question) Thanks, I have a question. That was a great talk by the way, that was good, thanks for the effort! I was wondering if you could, while supplying the command line arguments, change the genome build, since it looks like you hardcoded the hd19 into the code.
(Answer) The thing is, that you can provide so in the iGenomes configuration file. There are several genomes, one of those is the hg19 but there are more. There are not only human they're also from mice and so on. If you can check the key and this way all the files that you need the fasta file, the genome fasta file, these genome sizes that I told you, they are automatically rendered by the pipeline. In the case that you don't have them or in the case that you are running a genome that is not there, you can provide these parameters to the pipeline and these files and it will run. It's just for simplicity that I include this genome in the in the command.

[20:01](https://www.youtube.com/watch?v=59I-4wB1z4c&t=1201)
If there are no questions then i would like to thank Jose of course and the Chan Zuckerberg Initiative for funding of these talks. As usual this talk will be uploaded to youtube and if you have any questions later on you can always come to the slack channel of chipseq or for bytesize and ask the questions there. Thank you very much.

</details>
