---
title: "Bytesize: nf-core/hgtseq"
subtitle: Francesco Lescai - Department of Biology and Biotechnology, University of Pavia
type: talk
startDate: "2023-03-21"
startTime: "13:30+01:00"
endDate: "2023-03-21"
endTime: "14:00+01:00"
youtubeEmbed: https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69
embedAt: "hgtseq"
locations:
  - name: Online
    links:
      - https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69
      - https://doi.org/10.6084/m9.figshare.22317877.v1
---

This week, Francesco Lescai ([@lescai](https://github.com/lescai)) is presenting us the nf-core pipeline hgtseq. nf-core/hgtseq is a bioinformatics best-practice analysis pipeline built to investigate horizontal gene transfer from NGS data.
It uses metagenomic classification of paired-read alignments against a reference genome to identify the presence of non-host microbial sequences within read pairs, and to infer potential integration sites into the host genome.

<details markdown="1"><summary>Video transcription</summary>

:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69&t=1)
(host) Hello, everyone, and welcome to this week's bytesize talk. I'm very happy to have with me today Francesco Lescai from the University of Pavia at the Department of Biology and Biotechnology. He is very, very busy in nf-core, and among other things, he also worked with Sarek, but today he's going to talk about another pipeline, which is nf-core/hgtseq, and off to you.

[0:27](https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69&t=27)
Thank you, Franziska. Today I'm going to give you a bit of a background for this pipeline and the motivation that inspired us to initiate this project. I'm going to describe the pipeline components. I'll give you some usage indications on the performance of the pipeline, and then I'll describe a bit of future perspectives, which is our homework, basically. I'm going to start with the acknowledgments here, first to Simone Carpanzano, who's the lead author of this pipeline, but as you might imagine, he's heavily engaged in preparing the defense of his Bachelor of Science now, so he couldn't present today. Mariangela Santorsola, who's a key person in my lab, and she's also contributed to the publication that describes this pipeline. Then this is very important, I think, because the value of the nf-core community is the availability of all the modules that we also have used in our pipeline, so a very important acknowledgement here is to all the authors of the different modules that we have used and which actually make the added value of nf-core so important.

[1:43](https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69&t=103)
Starting from the background of this pipeline, horizontal gene transfer. This is a very known and studied process in biological organisms, and it refers to the transfer of genetic material between two different species when they are in close proximity. This has been very important in evolution because it has contributed to new traits, it creates adaptation to new environments, and also the capability to use new sources in different organisms. It's been crucial in the evolution, as I mentioned, particularly in archaea and bacteria, but not very much has been known about this phenomenon happening in higher organisms like mammals, for example.

[2:36](https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69&t=156)
Our motivation was mostly inspired by a paper several years ago that described the existence or the detection of microbial reads in exome sequencing data in human projects. That paper was really inspiring for us in the sense that it did highlight that microbial sequences have been found in exome sequencing data, which means the coding part of our genome, and it did open a huge lot of questions about these phenomenon in higher organisms, and it definitely needed end-to-end tools to investigate what is happening there. Of course, I put here a funny picture of the microbiome because if you remember the definition I just gave, which is transfer of genetic material between species that are in close proximity, then we and many other mammals are the leading example of this close proximity between different species, and we have a whole set of microorganisms that live with us and contribute to our own biology. Clearly, there's a lot to investigate here.

[3:53](https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69&t=233)
A couple of definitions for the pipeline that we have developed. First of all, when you map next generation sequencing reads to a host genome, and in our example, a human genome, you could have several scenarios. The first scenario, which is the most common, is that if you do paired end sequencing, both mates in the pair map correctly to the host genome. But you can also have a couple of additional scenario. One where only one of the two members of the pair maps and the other is unmapped, and one where both reads in the pair are not mapped to the host genome. We needed a definition for the pipeline, so we have identified these pairs where only one read is mapped to the genome and one is not, as "single unmapped". Then we have defined those where both members of the pair are unmapped as "both unmapped". You will find that these short definitions later on recurring in the picture and the slides that we present in a moment.

[5:05](https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69&t=305)
Of course, the importance of the pair where one mate is mapped and the other is unmapped is that it allows us to make assumptions about a potential integration site. Because of course, we can measure and evaluate the abundance of taxonomic IDs from every read that is not mapped to the host genome. But for those that are members of a pair where one of the two is actually mapped, we can additionally try to infer where that potential integration has happened thanks to the coordinates that we have from the mapping of the mapped member of the pair.

[5:48](https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69&t=348)
This is the pipeline overview. The pipeline I think is relatively straightforward and includes a part dedicated to the alignment of quality control, then the conversion and parsing of the reads that I just illustrated and classification using Kraken at the moment, and then a last phase of reporting. We're going to see each of these steps in the pipeline in a moment. The pre-processing is very important because it's been designed to be plugged downstream to other studies. I made the example of the initial paper that inspired us to develop this pipeline. That was the discovery of microbial reads within human exome studies. Our own idea, particularly because we have also contributed to Sarek, was to plug this type of pipeline downstream into those kinds of pipelines like Sarek. Accepting the bam files of the alignments that have been produced by human exome or whole genome sequencing studies, and then use the pipeline to process all those reads that have not been mapped. But the pipeline also starts from a fastq, so using raw reads, and it does a standard alignment to the host genome using BWA.

[7:10](https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69&t=430)
[interrupted stream]

[7:39](https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69&t=459)
We do this using samtools and using the bitwise flag, 13 and 5. Then we further parse the potential integration for the single unmapped reads using the information from the mapping coordinates of the mapped member of the pair. At this moment, we are using Kraken 2 to classify taxonomically the reads, and in particular we have chosen this tool because we're using the k-mer classification that is given as a sliding window in the NGS read that we are analyzing. Also as a way for interpreting the results and doing further QC on the outcome of the taxonomic classification.

[8:26](https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69&t=506)
These all goes into a reporting phase of the pipeline. We generate traditional Krona plots that are generated per group. If your analysis has one, two, or three different groups, we group the sample of Krona plots per category of your samples. We obviously use MultiQC for the reporting. This also includes the classification of a view of all the reads, thanks to the parsing of Kraken 2 outputs. Then we perform a preliminary analysis using RMarkdown with a parameterized RMarkdown file, which also adds a couple of important information to the preliminary analysis. One is a classification score. We try and use the information that Kraken 2 gives us in the output, in order to give a classification score to each of the reads to further allow us to filter based on the quality of the taxonomic classification. Important information here is how much of the read has been classified and has been assigned to that taxonomic organism, which appears in the result. Then we have also curated from a number of publications a list of contaminants that are known to affect DNA extraction and DNA extraction kits. We have further classified the contaminants depending on their potential role in human diseases as well, because we are particularly interested in analyzing these phenomena in humans.

[10:05](https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69&t=605)
A couple of indications about the usage, this is a typical common line to start the pipeline. We will use the input sample sheet as a comma separate value as most of the nf-core pipelines. Then we use the iGenomes genome indication, we use the host genome there, so this is the first part that performs the host genome alignment. Then we pass on the host taxonomic ID which is used to filter the results in the R Markdown report. Two very important parts of this common line are a path to the Kraken database and a path to the Krona database, which can be either indicated as a path if you have it locally, or as a tar.gz file, which can also be online or in a repository that you might have in a cloud resource. The inputs, as I mentioned in the beginning, can be either raw reads with a FASTQ input as you can see in the first example, or already pre-aligned BAM files that are coming from another pipeline that you can see in the second example of input. Here I also have to say that the database for Kraken is obviously crucial for the classification because the whole point of this pipeline is in assigning a taxonomic classification to the un-mapped reads. The way the Kraken database has been built obviously will have a huge effect on the results that you're able to report. On the taxonomic IDs that you're able to detect in your reads.

[11:56](https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69&t=716)
A couple of words about the performance. We have tested this pipeline on different species, to demonstrate the existence of the phenomenon in not only humans, but also in other mammals. This is an overview of the execution of the pipeline on 10 exomes from humans. You can see that they are executed in our local cluster in about three hours. This is quite good. The pipeline is very smooth in its run. Then we have also reported CPU and memory usage for the most intensive tasks. There's nothing major to discover here, I mean, in particular, in terms of memory. Kraken and QualiMap are also quite intense. Again, the amount of memory that is used by Kraken definitely depends on the database that is used for the classification. QualiMap is known to be quite greedy with the memory. In Sarek, it has been swapped with mosdepth. We might do the same in a future version of the pipeline for the same reasons.

[13:10](https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69&t=790)
Homework, mostly. As I mentioned, Kraken has a very useful type of output where you can appreciate the assignments to taxonomic IDs by a sliding window of the k-mers the reader has been splitted into. This will allow us to draw much more information in terms of classification filters or heat maps that will allow us to investigate better the biology, especially regulating this type of events. We will probably dedicate some work to the optimization of the computing part of the pipeline. I just mentioned the issues with QualiMap. Certainly improvement on the preliminary analysis report, which is currently running only on humans. Also consider the introduction of alternative taxonomic classifiers. Here we have a number of examples in other nf-core pipelines. I hope this is enough of an overview. For now, we have published a paper on the International Journal of Molecular Sciences very recently, where the nf-core community is a collective author in the publication as well. There you can find more details about, particularly about the scientific findings that we have collected by analyzing the different species we have used for testing of the pipeline. I'm open to take any questions.

[14:47](https://www.youtube.com/watch?v=nDaRt2L-tRw&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=69&t=887)
(host) Thank you very much. I have enabled now for anyone to unmute themselves if they want to ask any questions. Or you can also write questions in the chat. I will then read them out. It seems that it was very clear to everyone. If there are any further questions, you can always ask them in the nf-core Slack. Or you can ask directly, Francesco, I assume.

(speaker) Yes, definitely. Both on Slack and via email.

(host) Then I would like to thank you again, Francesco, and also the Chan Zuckerberg Initiative for funding our bytesize talks, and all of the audience for listening. Thank you very much.

(speaker) Thank you.

</details>
