---
title: 'Bytesize 22: nf-core/eager'
subtitle: James Fellows Yates - MPI-EVA, Germany
type: talk
startDate: '2021-10-05'
startTime: '13:00+02:00'
endDate: '2021-10-05'
endTime: '13:30+02:00'
embedAt: 'eager'
youtubeEmbed: https://youtu.be/fObuLeGhQRo
locationURL:
  - https://youtu.be/fObuLeGhQRo
  - https://www.bilibili.com/video/BV1q34y1U7Ak/
  - https://doi.org/10.6084/m9.figshare.16750387.v1
---

# nf-core/bytesize

Join us for a special pipeline-focussed episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 22: nf-core/eager

This week, James Fellows Yates ([@jfy133](https://github.com/jfy133/)) will tell us all about the nf-core/eager pipeline.

nf-core/eager is one of the largest and/or most complex nf-core pipelines. It is a best practise (meta)genomics primarily aimed at (but not limited to) processing ancient DNA. This talk will:

- Give an overview of the pipeline.
- Discuss the impact of the particularlities of ancient DNA.
- Cover some of the challenges and solutions encountered when writing this pipeline.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:39](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=39) So I’m going to talk about nf-core/eager, which is a best practice analysis pipeline for NGS sequence-based ancient DNA data analysis. During the talk today, I will briefly introduce what palaeogenomics is and why we need a special pipeline for it. I’ll also give a very brief overview of nf-core/eager and talk a little bit about some of the development challenges we had for this pipeline during DSL1.

[1:22](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=82) So what is palaeogenomics? It is a very diverse field that is focused on ancient DNA but there are many different facets to it. For example, there is a lot of work in human genomics, which is used to study past human history, there’s also animal - particularly megafauna palaeogenomics for studying past ecology and evolution, and then there is microbial genomics - particularly for pathogens, which allows us to study infectious diseases in the past. More recently, researchers including myself have started to study ancient microbiomes, and this also has several facets such as disease, human behaviour, and so on. Finally, there’s also a lot of work now on studying sediment DNA for ecology and evolution studies.

[2:18](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=138) Palaeogenomics is not particularly new and there have been previous pipelines dedicated to it. nf-core/eager is not a new pipeline; it is a new variant or a refactored version of a pipeline developed by [Alex Peltzer](https://github.com/apeltzer) about five years ago. This is the diagram that was in the publication, and as you can see here, there are three main sections: pre-processing, read-mapping, and genotyping and optional consensus calling.

[3:04](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=184) Now, you might be thinking the three stages are pretty standard; nf-core/sarek also does that. This difference here is the ancient DNA is very fragmented, damaged, and likely contaminated with modern DNA, and this complicates things in a variety of ways.

[3:34](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=214) To visualise this a bit, with fragmentation, I meant that over time with no repair mechanisms, DNA will be exposed to various chemical and other processes, which causes it to fragment. This makes them very short and in general down to about 30 bp, which is so short that you have very low specificity. So your reads can map to different regions in the genome or in the case of metagenomics, it can map to several different species. Because of this fragmentation, many can disappear over time. Libraries will have very few ancient DNA fragments resulting in low coverage and not very high confidence in variant calling and so on. On top of that, we have damage. Fragmentation does not always end in clean breaks on both strands, and sometimes this leaves you with this single-strand overhangs which exposes the nucleotides on the strands that still exist to other degradation processes such as deamidation that converts cytosines to uracils. So in addition to low coverage, we also have these substitutions that further complicate variant calling. Then in addition, there’s contamination where DNA in a buried environment comes into contact with other external sources of DNA such as other organisms in the environment, individuals handling the specimens, etc. Trying to distinguish between these, adds another layer of complexity.

[6:14](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=374) While all these things can make working with ancient DNA challenging, there is a range of authentication criteria that help separate true ancient DNA from the organism you want to study vs the rest. This includes things like damage profiles, fragment length distributions, edit distances, and the metagenomic component such as profiles in ancient human faeces etc.

[7:22](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=442) Palaeogenomics has been around since the start of the NGS revolution, so maybe 15-20 years. Since then, there have been a lot of developments in how to extract ancient DNA, so we can get pretty good samples. However, it can sometimes be too good, so particularly in the human genetics field, they are now routinely dealing with datasets of hundreds or thousands of samples per publication. The issue with previous pipelines like [EAGER [v1]](https://github.com/apeltzer/EAGER-GUI) was that they were not designed for HPCs, so they were very linear or were only able to work on a single node. On top of this, the field is becoming more and more interdisciplinary so there are a lot more diverse questions that are being asked such as whether ancestral populations died of a disease etc.

[8:43](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=523) The solution we came up with is converting the original eager pipeline and extending it into Nextflow, and we did this within nf-core/eager.

[8:59](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=539) I’m going to briefly go through this. If you have more detailed questions, get in touch with me on Slack. We have the basic genomics pipeline steps, and we’ve extended this with statistical analyses and a metagenomics screen for pathogen detection and authenticating reads.

[10:40](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=640) This is another representation of the graph showing just how complex the pipeline is. We have three different main routes: eukaryotic nuclear genome, microbial and smaller organellar genomes, and the metagenomic component where we do the taxonomic profiling. So it’s a very complex pipeline.

[11:33](https://www.youtube.com/watch?v=fObuLeGhQRo&list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&index=23) And this is why I want to talk a little bit about the main development challenges for the DSL1 version. I am hoping that some of these will be alleviated for DSL2, but that remains to be seen.

[11:45](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=705) So the first issue is complex input data; there are many different weird and wonderful ways in which people generate data that would go into the pipeline. So for example, it is possible to overcome sample damage in the lab through using different methods during library preparation. This can vary, and one would require different settings within eager to process that accordingly. Then there are several different sequencing configurations. Finally, there are heterogeneous input files.

[13:06](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=786) Our solution here was adopting the TSV input. We adapted heavily from sarek, but more importantly we made extensive use of something I call rerouting, which is basically using channel branching, filtering within that, and then merging. We had to do this in every single main step of the pipeline because of the variety of ways in which people provide their input into the pipeline. An example of the code of this adaptation can be seen on my screen for poly-g trimming, which is an artefact of NovaSeq or Nextseq data but does not apply to HiSeq data. This made the pipeline quite complicated, but was relatively quite efficient once we got the hang of it.

[14:25](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=865) As I said before, palaeogenomics can be a very interdisciplinary field. We can have people with several different specialisations wanting to play with data. Ancient DNA scientists are really good with sharing their raw data, so anyone can try this out. But what that meant is that we had the challenge of trying to communicate how to run the pipeline, train people to run the pipeline. The field of ancient DNA is not large enough to have courses on MOOCs or at universities etc.

[15:16](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=916) So the solution here was documentation! Since there was a lot of it, we thought it could sometimes be a little dry and not too easy to read. So the question was how we could make it more interesting and engaging to read and follow. So I collaborated with [Zandra Fagernäs](https://github.com/ZandraFagernas) ([Bytesize#14](https://nf-co.re/events/2021/bytesize-14-graphic-design) to draw these schematic images that included notes on how to interpret the outputs, especially since ancient DNA can be different as I mentioned earlier. In addition, we tried to construct it in a way in which it could be reused for training purposes. We hope that by doing this that a lot of students that begin working with ancient DNA and the nf-core/eager pipeline will be educated before they run it.

[17:18](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1038) The final thing was that there were a lot of opinions but no real standards in the field. Ancient DNA can be very competitive, particularly human population genetics and there are several big labs working in the field. But what this meant is that there were constantly changing standards and settings. This made it difficult for us to know which tools and parameters to use to analyse this kind of data.

[18:12](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1092) Our solutions to this started quite early on. [Alex Peltzer](https://github.com/apeltzer) added a little tool that helped estimate the optimal mismatch rate for [`bwa aln`](http://bio-bwa.sourceforge.net/bwa.shtml) during seeding. So developing a little tool ecosystem around this that feeds in and out of eager. To cope with the changing standards, I use Twitter, where I use the 'hive mind' to assess the general consensus that allows us to improve the various features in the pipeline. I’ve also repeatedly introduced it at workshops; this is very important to capture people from different fields in order to get their feedback on the pipeline.

[19:49](https://youtu.be/fObuLeGhQRo?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1189) So to summarise, palaeogenomics is complicated, but it is a fun challenge. Broad documentation helps in interdisciplinary fields. Finally, be active in outreach, not just in support, but to try and sell your pipeline because it keeps your project alive past publication and nf-core/eager will hopefully continue to be used and developed by scientists in the ancient DNA field.

</details>
