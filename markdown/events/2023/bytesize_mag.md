---
title: 'Bytesize: nf-core/mag'
subtitle: Sabrina Krakau - University Tübingen, QBiC
type: talk
start_date: '2023-02-28'
start_time: '13:00 CET'
end_date: '2023-02-28'
end_time: '13:30 CET'
youtube_embed: https://www.youtube.com/watch?v=IiorfDHeoLo
embed_at: 'mag'
location_url:
  - https://doi.org/10.6084/m9.figshare.22210879.v1
  - https://www.youtube.com/watch?v=IiorfDHeoLo
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: nf-core/mag

This week, Sabrina Krakau ([@skrakau](https://github.com/skrakau)) is going to introduce nf-core/mag. nf-core/mag is a bioinformatics best-practise analysis pipeline for assembly, binning and annotation of metagenomes.

<details markdown="1"><summary>Video transcription</summary>
**Note: The content has been edited for reader-friendliness**

[0:01](https://www.youtube.com/watch?v=IiorfDHeoLo&t=1)
Hello, everyone, and welcome to this week's bytesize talk. I'm happy to present to you today Sabrina Krakau. She is situated at QBiC at the University of Tübingen. She is talking today about the nf-core pipeline mag and off to you.

[0:21](https://www.youtube.com/watch?v=IiorfDHeoLo&t=21)
Thanks Franziska for this kind introduction. I'm very happy that it finally works out to also present the nf-core/mag pipeline to all of you. This pipeline you can use for metagenome hybrid assembly and binning. The goal of this pipeline is to analyze microbial communities by recovering individual genomes. This might be, for example, particularly useful if you do not have a complete set or high quality reference genomes given. Such microbial communities could be everything, for example, environmental samples, but also host associated communities such as the gut microbiome.

[1:02](https://www.youtube.com/watch?v=IiorfDHeoLo&t=62)
The microbiome samples can be processed with metagenome shotgun sequencing, which generates short reads. The nf-core/mag pipeline then essentially combines these reads and assembles them to larger contigs. In a fine downstream genome binning step, it bins these contigs to so-called metagenome assembled genomes, or also called MAGs. These MAGs can then further be annotated and also taxonomically classified. That's the concept of the nf-core/mag pipeline. As for many nf-core pipelines, the development of this was a quite large community effort with many different contributors, so just mentioning the main important ones. It was started by Hadrien Gourlé, then Daniel Straub contributed a lot since early on when I joined, and also since last year, James Yellows Yates is the main contributor of this pipeline.

[2:04](https://www.youtube.com/watch?v=IiorfDHeoLo&t=124)
Now I would like to mention the key features of this pipeline. It can perform a hybrid assembly using both short Illumina and long nanopore reads. This is useful because if you have assemblies generated only from short reads, they are often highly fragmented. By using additionally longer reads, this can improve the contiguity of such resulting assemblies. The pipeline also performs a genome binning step and optionally also a binning refinement step, then can taxonomically classify the resulting bins and also provides a comprehensive QC statistics. Furthermore, it can utilize sample-wise group information. This can be used for the core assembly. This is important if you have data sets where you know that certain strains are present across multiple samples, such as within longitudinal data sets. Because the core assembly can improve or increase the sequencing depth, this also allows to recover more lower abundant genomes. Additionally, the group information is also used for the computation of core abundances, which is used in the genome binning step. Furthermore, the pipeline also allows the handling of ancient DNA, because it's containing ancient DNA validation sub-workflow, which is rather specific for this pipeline. A previous version of this pipeline was already published at the beginning of this year in NAR Genomics and Bioinformatics, so if someone's interested in more details, you can also have a look at this application note.

[3:44](https://www.youtube.com/watch?v=IiorfDHeoLo&t=224)
Here you can see an overview of the pipeline. The pipeline starts with different pre-processing steps and QC, then the actual assembly is performed with a final genome binning step. Here in green you can see the processes or different tools that are run by default by this pipeline. In the following I would like to guide you through the different steps of this pipeline in more detail. Just first, how can we actually run it? So here you can see an example of the Nextflow command that is typically used and in order to run it with default settings, just provide a sample sheet as input file.

[4:26](https://www.youtube.com/watch?v=IiorfDHeoLo&t=266)
Here you can see an example how the sample sheet looks like for this pipeline: it contains five columns. The first column contains a sample name, the second column contains a group name, in this case all samples belong to the same group. Then you have to provide the path to the input read files, either only to the short read or to the short and long read, so the long reads are optional. Starting with this sample sheet file now, or if you have only short reads you can also just provide a fastq file directly. The pipeline then pre-processes the short and long reads separately from each other with different pre-processing steps. I do not want to discuss them in detail. Maybe just mention that the host reads can also be removed by mapping the reads to given reference sequences. This information is also used indirectly for the long reads, since the long reads are filtered based on the already filtered short reads. The short reads can then further be taxonomically classified already. This can serve for example as a quality control in order to check for potential contaminations.

[5:41](https://www.youtube.com/watch?v=IiorfDHeoLo&t=341)
After these pre-processing steps then the actual assembly is done. This can be done sample-wise or the group information can be used in order to run a whole assembly, however by default this is done for each sample individually. By default the tools SPAdes and MEGAHIT are run both. However, you should keep in mind that if you have long reads given and you are interested in the hybrid assembly then only the tool SPAdes can be used for this. Then the tool QUAST is used in order to assess the quality of the resulting assemblies and also the assemblies are further processed with the tool PRODIGAL which predicts protein coding genes for this.

[6:26](https://www.youtube.com/watch?v=IiorfDHeoLo&t=386)
That's the assembly part and the contigs of this assemblies are then further processed in the genome binning step, where the tools MetaBAT2 and MaxBin2 are used, which now bin the contigs to retrieve the actual genomes. The results of these tools can also additionally be combined in a binning refinement step, which makes use of DAS tool. The quality of this bin is as well assessed with the tool QUAST and in addition the tool BUSCO is used which makes use of single copy orthologs in order to estimate the contamination on the completeness of the retrieved genomes. Additionally the pipeline also uses a custom script, which estimates the abundance of the individual bins, because it's also a relatively important output of this pipeline. Further downstream processes then the bins are further taxonomically classified by default using the tool GTDP-Tk, and also annotated with the tool PROKKA. Finally a multiQC report is generated and also a relatively comprehensive MAG summary report.

[7:41](https://www.youtube.com/watch?v=IiorfDHeoLo&t=461)
How does the output of the pipeline look like? Besides all the individual results part of the individual tools, the pipeline generates a clustered heat map showing the MAG abundances across different samples. You can see an example how this looks like and if you would see here for example that certain samples clustered together, for which you know that they are originating from different groups. This might indicate that something has gone wrong. The pipeline also outputs the MAG summary, which I already mentioned. This contains for each bin for each MAG, the abundance information across different samples. The QC matrix from the BUSCO results and the QUAST results, and also taxonomic classifications from the tool GTDP-Tk.

[8:33](https://www.youtube.com/watch?v=IiorfDHeoLo&t=513)
And with this, I've shown you the rough overview of the pipeline and next I would like to show you some of the impact, different assembly settings can have. For this I simulated some mouse gut data set in the past with the tool CAMISIM, and I generated hybrid data containing Illumina data and Nanopore reads and generated two groups. Each with a time series of four samples. This might be the ideal case where a core assembly might be useful. Now I would like to show you some of the resulting assembly metrics that are commonly used.

[9:14](https://www.youtube.com/watch?v=IiorfDHeoLo&t=554)
Here you can see for example the total length of the resulting assemblies and then compared for different pipeline ones were different assembly settings were used. The lower two pipeline settings correspond to a sample-wise assembly and using either only short or short and long read, so hybrid data, and the upper two settings correspond to a core assembly. Again with short or short and long reads, and what we can see is that the total length of the resulting assemblies significantly increased both by using the hybrid setting, and by applying the core assembly setting. Similar results we also see when looking at the number of MAGs, so the number of genomes that could be retrieved from this data, and also when looking at the L50 values. This indicates that the actual setting that is used for the assembly within this pipeline can have a relatively huge impact on the results. It's definitely good that the pipeline provides different settings, so that you can really choose the correct setting for input data, and it might also be worth to compare different settings.

[10:32](https://www.youtube.com/watch?v=IiorfDHeoLo&t=632)
Another topic I would like to shortly mention is the resource requirements, because this came up quite often in the Slack channel, and it's also somehow difficult to estimate in advance, because it really differs depending on the input data. The main requirements are both for memory and time, coming from the assembly step. As I mentioned already it really differs for different input data sets and I collected some numbers just to give you a rough idea for different pipeline ones that were run by Daniel Straub on our compute cluster. For one rather small sample, which was a cultural sample, both MegaHIT and SPAdes required less than 25 gigabytes and were finished in a couple of hours. However, for a larger river sample data set, MegaHIT took already more than 100 gigabytes of RAM, and it took more than one day to finish, and SPAdes even took more than 900 gigabytes of memory, and it required more than nine days. There was another very large data set containing 15 soil samples for which also a core assembly was performed and for this MegaHIT required one terabyte and more than 17 days, and SPAdes could not even be run because it would have required more than two terabytes of memory.

[11:55](https://www.youtube.com/watch?v=IiorfDHeoLo&t=715)
This just shows that even for smaller data sets, you cannot run this on your laptop. In general, one can say that it depends on the sequencing depths, the number of samples, the complexity of the analyte metagenome, and also on the applied tool and setting. For this it might be worth noting that, that both assembly tools are run by default but MegaHIT requires much fewer resources than SPAdes, and if you do not want to compute a hybrid assembly it might make sense to consider the `--skip_spades` parameter. Additionally, the core assembly also increases the required resources because it pools samples. At least for one individual task, the required memory and time is much higher. This is something important to keep in mind, because also if you want to run it on larger data sets, you might want to provide a custom config file in order to adjust the resources required for your particular data set.

[12:53](https://www.youtube.com/watch?v=IiorfDHeoLo&t=773)
With this we have seen how we can run the nf-core pipeline for modern metagenomic data sets. As I mentioned already at the beginning, it can also handle ancient DNA. For this James and Maxime added an ancient DNA validation sub-workflow. This is particularly interesting because, as we know at least, there's no other such pipeline which can handle ancient DNA. What this essentially does is that it performs identification of possible ancient contigs by modeling ancient DNA damage patterns, and then polishes the contigs in order to remove the errors that are caused by the presence of such ancient DNA damages in order to allow more unpaired downstream analysis. This might be interesting for some of you to know that this pipeline can also handle ancient metagenomic data analysis.

[13:52](https://www.youtube.com/watch?v=IiorfDHeoLo&t=832)
With this, I'm already at the end of my presentation, just a few words on the outlook. The next release James already prepared, it just requires one more review. It contains another optional binning tool, namely CONCOCT. It will also allow optionally the bin QC with CheckM and GUNC. For the midterm future, it would be also very nice if a functional annotation step could be added, so depending on the strategy, either using HUMAnN 3 or eggNOG, and also a standalone long read assembly option would be very nice by using, for example, the tool meta-flye, such that the pipeline could be also run without short read data.

[14:40](https://www.youtube.com/watch?v=IiorfDHeoLo&t=880)
In general, if you are interested in contributing, or if you have any questions or problems you would like to discuss, you can join us in the nf-core Slack channel dedicated to the MAG pipeline, or have a look at our GitHub repository. We're always happy about feedback or particular bug reports and issues. With this, I would like to thank you for your attention. Then, in particular, my colleagues from QBiC, importantly Daniel Straub for many contributions, James and Maxime from the MPI for Evolutionary Anthropology, Hadrien, of course, and importantly, the whole nf-core core team and community for helping with the development, for reviewing, testing and creating issues. With this, I'm happy to take any questions.

[15:30](https://www.youtube.com/watch?v=IiorfDHeoLo&t=930)
(host) Thank you very much. There is indeed one question already in the chat.

(question) It was at the very beginning when you were talking about examples, and you mentioned CAMISIM. Could you explain more in detail what this is?

(answer) This is a tool which was also used in the CAMI challenge to simulate metagenomics data. It's using as input different genome sources. I used in this case a set of mouse genome sources, which was given from some mouse gut data sets. Then it can simulate Illumina and nanopore data and simulate also different taxonomic profiles. But the more details, I would also have to look up, it was quite a while ago. Was there any particular question about this?

(question cont.) No, it was just a question, "what is CAMISIM?", but I think James has now added some links to articles. If anyone is interested, they can have a look at that.

[16:43](https://www.youtube.com/watch?v=IiorfDHeoLo&t=1003)
(host) For anyone else, if there are more questions, you can now unmute yourself and just ask them straight away. Or you can put them in the chat and I will read them out for you.

(question) I would actually have a question. What happens to multi-mappers? I can imagine that if you have related bacteria that it would also map to different ones. How does the pipeline deal with that?

(speaker) I mean, this is handled by the assembly tools then somehow.

(question cont.) But are they removed or added to all of them? Any idea?

(speaker) Someone of the others are more in the details of this algorithmic parts of the assembler.

(audience) Do you mean when you're mapping back to the contigs or during the assembly itself?

(question cont.) During the assembly. I mean, you map to the genomes, I guess?

(audience cont.) No. We need to explain the main concept there. But there's some fancy maths magic that goes on which estimates which reads most likely go with each other based on the number of mutations they have with each other. There's some weird maths stuff which works out which is the best grouping.

(question cont.) Okay, then I misunderstood that part. Thank you.

[18:15](https://www.youtube.com/watch?v=IiorfDHeoLo&t=1095)
(host) Are there any more questions from the audience? It doesn't seem so. If you have any more questions later on, as you mentioned, you can always go to nf-core Slack and ask questions there. If this is now all the questions answered so far, I would like to thank Sabrina again for this very nice talk. Of course, as usual, I would also like to thank the Chan Zuckerberg Initiative for funding the bytesize talks and of course everyone in the audience for listening. Thank you very much.

(speaker)) Thanks.

</details>
