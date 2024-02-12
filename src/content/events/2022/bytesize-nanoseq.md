---
title: 'Bytesize: nf-core/nanoseq'
subtitle: Yuk Kei Wan - Genome Institute of Singapore and National University of Singapore
type: talk
startDate: '2022-06-21'
startTime: '13:00+02:00'
endDate: '2022-06-21'
endTime: '13:30+02:00'
embedAt: 'nanoseq'
youtubeEmbed: https://www.youtube.com/watch?v=KM1A0_GD2vQ
locationURL:
  - https://www.youtube.com/watch?v=KM1A0_GD2vQ
  - https://doi.org/10.6084/m9.figshare.20115506.v1
---

This week, Yuk Kei ([@yuukiiwa](https://github.com/yuukiiwa)) will talk about the newest developments in the nf-core/nanoseq pipeline.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=1)
Hello, everyone. I'm Franziska Bonath. I'm today's host, and with me is Yuk Kei. She is giving an introduction to nf-core nanoseq. Welcome.

[0:15](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=15)
Thanks for the introduction. I'm Yuk Kei, and I've been maintaining nanoseq for the past two years since I picked it up from Harshil Patel, Laura Warden, and Chelsea Sawyer, and also training. So I've seen it through its DSL1 days, converted it to the initial DSL2 syntax, and updated it to the newest DSL2 syntax with a lot of help from Chris Hakkaart. Chris is from Seqera lab, and he was previously in the University of Tübingen. If you have any questions on nanoseq, you can always reach out to us.

[1:05](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=65)
So without further ado, I'll just tell you a little bit about nanoseq. Nf-core nanoseq is a bioinformatics analysis pipeline for nanopore DNA and RNA sequencing data that can be used to perform base calling, demultiplexing, QC alignment, and also downstream analysis. In this bytesize talk, I will briefly introduce you to what is nanopore sequencing and why we need a specific pipeline for nanopore sequencing data. Through Slack channel conversations, I pretty much have found people having trouble understanding how to run nanoseq itself. So I'll just go through the basics for how to run different parts of nanoseq and also talk to you about some latest additions to the pipeline itself, which we are very excited to introduce you to.

[2:08](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=128)
I don't know how familiar the audience is with nanopore sequencing. Nanopore sequencing is a sequencing technology that's provided by Oxford Nanopore Technologies. There's this string of nucleic acid going through a pore, and this is called the nanopore. As it goes through a pore, current signals are emitted from the nanopore itself. These current signals can be translated into their specific nucleotide bases.

[2:46](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=166)
Nanopore sequencing most people know as one of the third-generation sequencing technologies, and it outputs long-read. The longest read is around 2.3 megabases. Nanopore sequencing was used in the telomere-to-telomere consortium, which completed the human genome finally. It is also used in identifying RNA isoforms because it can sequence full-length RNAs.

[3:22](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=202)
Nanopore sequencing is different because it has these kind of current signals that is outputted by the pore. These kind of current signals are not available in other sequencing technologies. With these current signals, one can use machine learning algorithms to extract biological information, such as DNA modifications, RNA modifications, poly-A tail links, and also RNA secondary structure, without having to do extra lab assays.

[4:01](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=241)
Here's nanoseq. It is relatively convoluted. Chris Hakkaart created this figure, which is so nice. It is color-coded, based on the kind of sample that you have. We have three different subway lines here. For the blue line, there's this DNA sample, and for the green line, you get the direct RNA, and also it's cDNA that's aligned to the genome. The orange line, is the direct RNA that's aligned to the transcriptome. In the subsequent slides, I will just make it bytesized. I will just talk about what kind of stuff we can do with these three lines.

[4:55](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=295)
For the first part of the pipeline, it involves base-calling, demultiplexing, QC, and also alignment. Base-calling starts from a fast5 directory containing a bunch of fast5 files. For the pipeline input you have to input it with the `--input_path` flag, with the fast5 directory. To correctly base call your sample, you have to specify the flow cell. Demultiplexing can start from either the fast5 directory or the demultiplexed fastq, where if you have a fast5 or a demultiplex fast5, you can demultiplex it with Guppy, and you can output demultiplexed fast5 files with ONT fast5 API. If you have a demultiplexed fastq file, you can demultiplex it with QCAT.

[6:13](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=373)
You also need to specify the barcode kit in order for it to be demultiplexed correctly. After the demultiplexing of the fast5, we implemented PyCoQC and nanoplot for quality checking the fast5 files. To demultiplex fastq, we have FastQC and also nanoplot for quality checking the fastq files. Alignment can either take in the fastq files from upstream processes or from the user input where the fastq file is already demultiplexed.

[7:06](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=426)
So for the blue line, which is the DNA, Chris added these DNA variant calling tools for DNA small variant calling and also for structural variant calling. You can choose between medaka and deepvariant for small variant caller. You can choose sniffles and cuteSV for a structural variant caller. The default is medaka for a small variant caller and sniffles for a structural variant caller. If you have more questions, if you have any questions on DNA structural variant calling, you can reach out to Chris. He knows a lot more than I do on this. He is also on the call, so he'll take questions if you are interested in these.

[8:07](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=487)
For RNA, if you have a cDNA sample or direct RNA sample, and if you align it to the genome, you can do transcript discovery and also quantification. These processes take in a sorted BAM from samtools. The default is BAMBU. BAMBU does both transcript discovery and also quantification. We also have another option, which uses stringtie2 for transcript discovery and featureCounts for quantification. After transcript discovery and quantification, if you have more than one group of samples, you can also do a differential expression analysis on the gene level with DESeq2 and on the transcript level with DEXseq.

[9:05](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=545)
On the green line, and this is also a new functionality we implemented, we included JAFFAL for detection of RNA fusions. It takes in a fastq file from either the sample sheet, where given that it is demultiplexed, and we can take it from upstream processes. You can start from fast5 files or demultiplexed fastq files.

[9:46](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=586)
The last part of the newest edition of the pipeline is the RNA modification detection. This one is a little different in a sense that per sample, you should have a larger directory. Within the directory itself, you should have fast5 subdirectories and a fastq subdirectory. Within the fast5 subdirectories, you need to include all the fast5 files. In the fastq subdirectory, please only include one base called fastq file. It goes through the alignment to the transcriptome, converts it to BAM, then prior to RNA modification detection, nanopolish is run for segmentation, and if you only have a single sample, you can detect M6A with m6anet, if you have multiple groups of samples, if you want to see the differential modification across the samples, you can run xpore, it does the differential modification analysis.

[11:13](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=673)
Just to summarize, here is the nanoseq pipeline. There are a lot of tools included, but you don't have to install anything other than nextflow and Docker, Singularity or Conda, depending on whether you are using an AWS cloud or you're using an HPC. With the latest release of nanoseq, it supports DNA variant calling, transcript discovery and quantification, RNA fusion detection, and also RNA modification detection. Thanks for listening, and I'm happy to take any questions.

[11:59](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=719)
(host) Thank you very much. Am I visible? I have to remove the spotlight. Anyway, I have now allowed everyone to unmute themselves if they want to ask questions. Otherwise, you can also put questions in the chat. There's actually a comment that we have in the chat from Olaitan, sorry if I butchered the name.

(question) He says that t2t completed _a_ human genome and not _the_ human genome. Also he thinks that sniffles2 exists, which is an enhanced caller for structural variants. Have you thought about sniffles2? I know this was not your main part, but...

(Chris, answer) I can probably jump in. I've seen that come out very recently, maybe in the last four months. So sniffles was initially added about 12 months ago. It's the first caller that we were interested in. In my opinion, it's also superseded by cuteSV, which when we did testing, it was actually the best caller. But in saying that, I haven't actually tested sniffles2 with a full data set, so I can't be sure if it's better or worse, but something we can definitely look into for a quick add in the future.
(speaker) Yeah. Thanks, Chris.

[13:34](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=814)
(question) Can you hear me? Hi. Thanks. Great talk. My interest lies in mostly native RNA transcriptomics, I wondered if there is a way to add poly-A tail measurements.

(speaker) Yeah, that's on our radar. We actually talk about adding the poly-A tail length detection to it, because one of our lab members actually added the poly-A tail functionality from nanopolish, but I'm looking into TailFinder right now. So yeah. Do you have any specific poly-A tail length prediction tools that you want to add in?

(question cont.) Well, nanopolish2 works really great with the native RNA, and TailFinder can do cDNA, however, it uses CPU, and it's kind of slower if you have a big data set, while nanopolish is much faster to do it, and I think there is a Shiny app where you can visualize it. And one more other thing is while we're at the three prime end, are there any tools available to add to this pathways for alternative polyadenylation characterization of some transcripts?

(speaker) Do you have any tools that you have specifically in mind on poly-A alternative polyadenylation? For long read?

(question cont.) For longer read at the moment, I am just looking at something called LAPA. It's on GitHub. I think the group is still working on a paper, but it's on GitHub. It's long read alternative polyadenylation and they are calling it LAPA.

(speaker) LAPA?

(question cont.) Yeah.

(speaker) I am aware that there are short read polyadenylation tools such as QAPA, LAPA, there are quite a few tools out there that do that. So I'm not exactly sure how translatable those tools are to nanopore reads, and so it will be great if you can suggest several long read tools that we can look into too.

(question cont.) I think there is a long, roundabout way of getting it through FLAIR (Full-Length Alternative Isoform analysis of RNA), and then there is, I think, Tapas has got something, but I think Tapas is more Pacbio orientated.

(speaker) Tapas as of like, because I'm aware of the TPPAS. Because I'm aware of the fact that there's another Tapas, which is TAPAS, all TAPs.

(question cont.) Yeah. No, it's not that one, the thing is this one is I think they've got their own kind of Java user interface, but underneath it, there is a lot of FLAIR and scanT3, that's where it does some of the transcript variants and poly-A detection and finding of A-B-A sites. But it'll complete the pathway quite nicely to take it from everything that you have currently to A-A-P-A and also poly-A tail, I think. Thank you very much.

(speaker) Yes, for sure, thank you for any suggestions.

[17:23](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=1043)
(comment) Yeah, I think I can just as a follow up to that, that the group that did the TPPAS, they did a very robust, would I call it software now, because I saw the demo by the PI of that group. I mean, it does a lot of things. If you're doing anything isoform related, you just need to, they've done so much work. Since you're doing transcript discovery, you can just figure out a way to maybe find some of the features that they have implemented in things TPPAS or there's another one, ORCAS, where you can integrate into your pipeline since you're already doing transcriptomics and structural variant stuff, you know, with long reads. So that would really help your pipeline to become more robust.

(speaker) Okay, great, thank you, thank you for these suggestions.

[18:18](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=1098)
(question) There's also another question in the chat, Ido is asking, do you plan to include genome assembly tools to the pipeline?

(answer) We actually talked about this during the lab meeting last week, and we hope to include Raven into the genome assembly, because I suppose it is also nanopore based. We are looking to include that, and do you have any specific tools that you have in mind that you would suggest us to include it into?

(host) Ido, you could also unmute yourself if you wanted to, otherwise you can write in the chat.

(question cont.) Yeah, thanks. Can you hear me? I don't have anything in particular in mind, but you know, just the common pipeline tools, you know, that whether it is in raccoon or nanopolish or, you know. It depends obviously on which organism we 're looking at. I would say the common pipelines that would either take nanopore or reads by itself, or whether taking hybrid, both Illumina and nanopore to do the assemblies. We're just looking at some of the common pipelines. But to put them inside this workflow will make it very easy to use unless there's already a pipeline that exists that is designated for genome assembly.

(answer cont.) Okay, great. Awesome. Yeah, thanks for these suggestions.

[20:16](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=1216)
(question) There is another question in the chat. What is the best practice process for nanopore sequencing metagenomics?

(comment) That's a complicated question, actually. What's your answer?

(speaker) Yeah, that's a good question.

(comment) I think I can jump in. Somebody did a benchmark of metagenomics tools, and really came up with a conclusion that there is no best. So I think that's the simple answer to that question, because some of the well-established ones didn't even perform well when he did this benchmark, and he presented these at a conference, I think about a month ago, that I was part of, it's a long read conference where they were just doing different kinds of presenting tools and stuff like that. In a way it was, when it came to the metagenomics part, when he did his presentation, there was no best metagenomics tool for long read sequencing.

(speaker) Yeah, interesting.

[21:29](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=1289)
(question) I have a question for the nf-core core team. Is there a pipeline for metagenomics? I suppose MAC is for metagenomics, right?

(host) Is someone from the core team here?

(answer) I was just checking that. I think there is something that's pretty close. I think in terms of nanoseq, doing metagenomics might be slightly outside the scope. I think one thing we've found developing this pipeline is that it's a bit of a beast already. You know, you could easily split this pipeline into three different pipelines, one's for DNA, one's for standard RNA-Seq, and one for these isoform detection. Yes, we could look at trying to include it, but I think it, personally, I think it might be a step too far, and I am suspicious that there is another pipeline that does some form of metagenomics, but I can't remember the name of it, so please don't quote me on that.

[22:36](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=1356)
(question) Hi guys, just one more comment or question about including short reads. I think that's what the TPPAS pathway includes, a lot of short reads for transcript variants. It could be beneficial for both the assembly, but also for more complete transcript variants, if there is an option to include Illumina short reads, and obviously it increases some depth if you want to do the DESeq differential expression analysis, because I don't know if the depth is high enough in some of the sequencing. So that's another comment. One more comment is about de-novo modified-based detection with Tombo. Can that be added? I know Tombo is not really... well.

(speaker) Yeah, we are exploring that right now.

(question cont.) Okay. Thank you very much, great work, again, quite excited about this pathway.

[23:43](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=1423)
(comment) There's a comment in the chat that talks about the scale or the scope of this work you're currently doing. He was just saying that, do you think that this workflow is becoming too big, and I was actually going to say the same thing. So what I was going to say is that you can actually focus for now, maybe, you don't have to take my suggestion, on maybe the RNA-based analysis and just make that as robust as you can, or focus on the DNA-based analysis, and then when you think that one is decent enough in terms of scale, you can then move on to the other one, you know, I think somebody is just making a similar comment.

(speaker) I think that's a good suggestion. I think what Chris and I are a good team, he does the DNA part of it and people from the lab I work at, we do the RNA part of it. I think we also have some concerns about when is the pipeline out of scope, and so definitely that's on our radar to think about. Chris, do you have anything to add?

(Chris) I completely agree. I think initially last year when we started adding all of these new features it made sense because the front end of the pipeline was more or less the same and it made sense to recycle it. But now with a lot of the pipelines building these really awesome subworkflows, which can be shared and integrated multiple pipelines, it would be nice to lean into that side of the community and share what we're doing and be shared with as well. One of the things about nanoseq is that the sample sheet is a little bit atypical in that you can specify genomes for different samples or different genome for different samples and different alignments and things like this. So nanoseq is already coloring a little bit outside the lines of your typical nf-core pipeline. It is something we've spoken about and I think we'll speak about it again very soon about trying to bring it back to that nf-core way of doing things. As a part of that, I could see that we may consider splitting the pipelines. But that's something we'll have to talk about soon, I think.

(speaker) For sure.

(comment cont.) Thank you.

[26:20](https://www.youtube.com/watch?v=KM1A0_GD2vQ&t=1580)
(host) Okay, if there are no more questions at this moment, I thank you again, Yu Kei. I also would to take the chance to thank the Chan Zuckerberg Initiative for funding these talks. If there are any more questions to anyone here, you can always come to the Slack channel for bytesize talks or specifically for nanoseq and ask your questions there and you might get an even more detailed answer. So thank you very much, everyone.

(speaker) Thank you.

</details>
