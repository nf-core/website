---
title: 'Bytesize: nf-core/taxprofiler'
subtitle: Sofia Stamouli, Karolinska Institutet
type: talk
start_date: '2023-01-17'
start_time: '13:00 CET'
end_date: '2023-01-17'
end_time: '13:30 CET'
embed_at: 'taxprofiler'
youtube_embed: https://www.youtube.com/watch?v=p1EQtidJiUY
location_url:
  - https://www.youtube.com/watch?v=p1EQtidJiUY
  - https://doi.org/10.6084/m9.figshare.21916416.v1 (video)
  - https://doi.org/10.6084/m9.figshare.21916386.v1 (slides)
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: nf-core/taxprofiler

This week Sofia Stamouli ([@sofstam](https://github.com/sofstam)) will talk about nf-core/taxprofiler, a bioinformatics best-practice analysis pipeline for taxonomic profiling of shotgun metagenomic data. It allows for in-parallel profiling with multiple profiling tools against multiple databases, produces standardised output tables.

<details markdown="1"><summary>Video transcription</summary>
**Note: The content has been edited for reader-friendliness**

[0:01](https://www.youtube.com/watch?v=p1EQtidJiUY&t=1)
(host) Hello, everyone, and welcome to the first bytesize talk of 2023, and I'm very, very happy to have Sofia Stamouli present today a new pipeline called nf-core/taxprofiler, which is soon to be released, I've heard. Off to you, Sofia.

[0:23](https://www.youtube.com/watch?v=p1EQtidJiUY&t=23)
Hello, everyone. I'm going to talk about nf-core/taxprofiler, which is using the GitHub description "a bioinformatics best practice analysis pipeline for taxonomic classification and profiling of shotgun metagenomic data". In the talk today, I will briefly introduce what is shotgun metagenomics and how the development of tax-profiler started. I will give an overview of the nf-core/taxprofiler pipeline and how you can use and run the pipeline, as well as our upcoming development plans.

[1:07](https://www.youtube.com/watch?v=p1EQtidJiUY&t=67)
To start with, what is shotgun metagenomics sequencing? I borrowed the description from Kim's paper from 2017 that describes shotgun metagenomic sequencing as the untargeted sequencing of all microbial genomes present in a sample. It allows for the determination of the taxonomic diversity in a sample. We may be looking at bacteria, viruses, fungi, archaea, or a combination of those, that are present in a sample. The development started in February 2022 by James Fellows Yates and Maurice Beber. We at Karolinska Institutet joined during the online hackathon in March.

[2:00](https://www.youtube.com/watch?v=p1EQtidJiUY&t=120)
With that, I would like to mention that this is really a community-based development. There are a few nf-core pipelines, like eager and mag, that support some sort of taxonomic classification. But they only support one classifier. Each classifier is tailored for specific purposes, each one has its own custom output format. There was really a need to have a pipeline that would support taxonomic classification and profiling of metagenomic reads using both, multiple tools and multiple databases. There are at the moment a few examples of how you can use nf-core/taxprofiler. Some of those different contexts is pathogen detection in clinical metagenomics. One can use it for a comparative microbiome diversity analysis as well as detection of food DNA from enzyme microbiome samples. But of course, they are not only limited to those.

[3:15](https://www.youtube.com/watch?v=p1EQtidJiUY&t=195)
This is the overview of how the pipeline looks like. I will go into more details in the next slides. To start with, it supports both short reads and long reads. The first step is the sequencing quality control. Right now, FastQC is used as a default. But during hackathon in October in Barcelona, falco has been added as a drop-in replacement, which supposedly is an improvement especially for long reads. The user can choose between either FastQC or falco. Next we have the pre-processing steps. All of those are optional and up to the needs of the user. We have dedicated tools for each sequencing technology. The first step is the adapter removal where fastp and AdapterRemoval is supported for short reads and Porechop for long reads. Then taxprofiler allows for removal of low complexity reads with BBDuk and PRINSEQ++ for short reads and Filtlong for long reads. The user can also choose to remove the host reads using bowtie2 aligner for short reads and minimap2 for long reads. As the last step of those pre-processing steps, taxprofiler allows for concatenation of multiple FastQ runs or libraries of a sample.

[5:08](https://www.youtube.com/watch?v=p1EQtidJiUY&t=308)
The last step of taxprofiler is, of course, taxonomic classification. Right now we support nine classifiers/profilers with kraken2 being paired with Bracken, KrakenUnique, MetaPhlAn3, MALT, DIAMOND, Centrifuge, Kaiju, and mOTUs. Each profiler can be executed with multiple databases. It's with their own settings. Each profiler has its own output. Because each profile classifier has its own output format, taxprofiler supports standardized and aggregated taxon count tables with the help of taxpasta, that is a Python package and with Moritz Beber is leading the development. It stands for taxonomic profile aggregation and standardization. I added the link to the GitHub repository.

[6:23](https://www.youtube.com/watch?v=p1EQtidJiUY&t=)
In this slide, I'm going to talk about how taxpasta works. Here you can see an example of how the output of the kraken2 classifier looks like. It has six columns: the percentage of reads covered, the number of reads covered, the number of reads assigned; This column here describes the taxonomic level, this one describes NCBI's taxonomy ID, and this is the scientific name of each taxon. This is how the output from the Kaiju classifier looks like. It has five columns, it also has header and it is very different from kraken2. This is the case for all the different classifiers. With taxpasta, we are really able to have a standardized output format for each classifier. The output format looks like this. It has two columns. The first one describes the taxonomy ID, and this column describes the read counts.

[7:38](https://www.youtube.com/watch?v=p1EQtidJiUY&t=)
About how to run the pipeline, one would need two input sample sheets: one describing the FASTQ files and one describing the databases. This is how format of the sample sheet that describes the FASTQ files should looks like. The first column should describe a unique sample name. The user can add a run accession, and should describe the name of the sequencing platform, as well as the path to the FASTQ files. Regarding the sample sheet describing the databases, this is how it looks like. It is four columns. In the first column one should give the name of the classification tool. Here is a unique name, based on the database. In this column, the user can specify the parameters that they would like to use. The fourth column describes the path to each database. About `TOOL1` and `TOOL2` (the argument here), those can be replaced by its classifier or profiler that is desired by the user. The last argument, the `perform_step`, this can be replaced by pre-processing or post-processing steps.

[9:26](https://www.youtube.com/watch?v=p1EQtidJiUY&t=566)
About our future plans, we would like to support more taxonomic classifiers, particularly for long reads. We would like to add an assignment validation step by aligning matched reads to identify the genomes, and we would like to add the workflow for database construction. But before we go on with the implementation of those plans, please stay tuned for the first release in January. With that, I would like to thank James Fellows Yates in Germany and Moritz Beber in Denmark, as well as my colleagues here in Sweden: Tanja Normark, Mahwash Jamy, Lauri Mesilaakso, and of course all the collaborators that contributed with different classifiers and issues in taxprofiler. If you have any questions, please reach out to our Slack channel with the hashtag taxprofiler, and that's it. I'm happy to answer any questions.

[10:43](https://www.youtube.com/watch?v=p1EQtidJiUY&t=643)
(host) Thank you very much, Sofia. Are there now any questions in the audience? You can either write your questions in the chat, or you can unmute yourself. I allowed that now for anyone. If there are no questions at the moment, I actually have a question.

(question) I was wondering why there are so many of these profilers, because, I mean, if there was one that actually would work properly, then you would only need that one.

(answer) The metagenomics field is very broad, and with those classifiers, they're based on different algorithms, and they cover different needs.

(question cont.) The final output that you have now, is that an average of what the different ones detect, or?

(answer cont.) We have a different output for each classifier, and we have, with the help of taxpasta, we are able to have a standardized output for each of those classifiers.

(question cont.) Okay, but you will get a separate output for each classifier?

(answer cont.) Yes. At the moment, yeah.

[12:08](https://www.youtube.com/watch?v=p1EQtidJiUY&t=728)
(question) Then we have here questions in the chat. One is from Juan. Do you have to download the databases manually?

(answer) Yes. We do not support it right now. It's in our future plans, maybe to add a workflow for database construction, but the user has to do it by themselves right now.

[12:29](https://www.youtube.com/watch?v=p1EQtidJiUY&t=749)
(comment) Then a comment from James. I guess it is for the profiler question I had. He says it's also a fun problem for computer scientists. Thank you.

[12:42](https://www.youtube.com/watch?v=p1EQtidJiUY&t=762)
(host) Are there any more questions? It doesn't seem to be like. If there are questions later on, you can always reach out, as you mentioned, in the Slack channel for taxprofiler, or also in the bytesize channel. Otherwise, I would like to thank Sofia again for this great talk, and of course, also, the Chan Zuckerberg Initiative for funding these talks. Thank you very much, everyone, and I hope to see you next week.

</details>
