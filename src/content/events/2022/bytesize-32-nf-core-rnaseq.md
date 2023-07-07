---
title: 'Bytesize 32: nf-core/rnaseq'
subtitle: Harshil Patel - Seqera Labs, Spain/UK
type: talk
start_date: '2022-02-08'
start_time: '13:00 CET'
end_date: '2022-02-08'
end_time: '13:30 CET'
embed_at: 'rnaseq'
youtube_embed: https://www.youtube.com/watch?v=qMuUt8oVhHw
location_url:
  - https://www.youtube.com/watch?v=qMuUt8oVhHw
  - https://doi.org/10.6084/m9.figshare.19161176.v1
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 32: nf-core/rnaseq

This week, Harshil Patel ([@drpatelh](https://github.com/drpatelh/)) will tell us all about the nf-core/rnaseq pipeline.

nf-core/rnaseq is a bioinformatics analysis pipeline used for RNA sequencing data.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=1)

(host) Hi, everyone. Thank you for joining in for today's bytesize talk. I would like to begin by thanking our funders and the Chan Zuckerberg Initiative for supporting all events. Just some preliminary information. This talk is being recorded and the video will be uploaded on YouTube and shared on Slack and our website. The talk will be about 15 minutes after which we will have a Q&A session where you are free to send your question in the chat box where it will be picked up from there or unmute yourself and ask your question. Today, we'll be having Harshil Patel, the head of scientific development at Seqera Labs, who will be presenting to us about the nf-core RNA-Seq pipeline, which is a bioinformatics pipeline used to analyze RNA sequencing data obtained from organisms with a reference genome and annotation. Over to you, Harshil.

[0:54](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=54)
Thanks Simeon. Good afternoon, everyone, and thank you for joining what is the 32nd bytesize talk of this awesome series. I'm Harshil Patel. I am head of scientific development at Seqera Labs. I'm also one of the long-term contributors to nf-core and various other pipelines that we have on nf-core. If you want to know more about me, there's a link here. Just click on that. It's a blog I wrote recently when I joined Seqera Labs.

[1:25](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=85)
Jumping directly into some numbers. This pipeline is one of the oldest and most popular pipelines on nf-core. The numbers are just staggering, and they always surprise me when I see them. We've got 400 forks, almost 60 contributors. It's also almost 700 people on Slack, and it's also one of the most active channels on Slack, where people are reaching out for help and coming to join to ask questions, and also just as a forum to discuss the pipeline as well. Over the years, this has really been one of the main pipelines that we've had on nf-core, and I would say that a lot of this has really been possible as a result of the testament to Nextflow itself, which is the underlying language that we're using. It's just allowed us to have access to communities, infrastructures, and other stuff that we wouldn't normally be able to do with a pipeline like this.

[2:38](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=158)
The pipeline itself has gone through various releases now over the years. As I mentioned, Phil from NGI initially pushed this when nf-core was first starting up and it was one of the main pipelines that he pushed here, and then it went through various iterations of updates, and Alex Pelzer got involved in between for a while, and then there was a sort of a gap for about a year where we really needed someone to sit down and update the pipeline. That's where I got involved, mainly in helping out with the implementation of the pipeline. Before, up to version 1.4.2, the pipeline was written in Nextflow DSL1, and then some of you may know that Nextflow has now a new DSL2, it's a more modular language. For us, I think that was the perfect opportunity to start from scratch, rewrite this pipeline essentially from scratch in DSL2 to allow us to have a proof of concept as to how it would work on nf-core, because obviously we want other pipelines to adopt similar syntaxes and principles.

[3:49](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=219)
I went about coming up with the first iteration of DSL2 at that point, and we released version 2.0. Since then, we've now changed and adopted the way that we're using DSL2, partly due to updates. Mahesh helped out with, what is the second iteration of DSL2 that we've now got on nf-core. It's constantly improving, it's being adopted more and more across nf-core, and you'll be able to see that in version 3.5. this pipeline has really become the cutting edge or the gold standard in terms of what we're doing with Nextflow implementations.

[4:31](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=271)
In terms of the RNA-Seq itself, it's probably one of the most popular applications of next-generation sequencing, and most people doing experiments will have come across some sort of RNA-Seq data, I imagine, especially bioinformaticians. What you're doing is you are quantifying the expression of genes in a genome at a given time. This is typical of bulk RNA-Seq sequencing. You then want to get a quantification of what your genes are, what the expression of your genes are like in one condition compared to another, and then figure out what is different and try and put that in some sort of functional context, like looking at pathways or doing further experiments to figure out how expression is impacting functionally what you are doing or how you're perturbing the cells.

[5:27](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=327)
A typical pipeline for this would be you have your reads, you do some cleaning of these reads by removing adapters and stuff that you get off the sequencing technologies, do some sort of QC. In this case, we don't actually have this bit in the pipeline, but it's probably something we may add later. I'm still thinking about how to do this properly. But this bit here allows you to sample reads and essentially automatically infers strand specificity and then plug that directly into alignment algorithms, which need this information. You would do some cleaning, and then you would map to the transcriptome, and then you can get some QC out from your genome BAM files as well, like looking for intronic rates or genomic contamination and all sorts of other really useful QC metrics from your genome alignments.

[6:19](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=379)
But most importantly, you also get the gene counts out, and this is essentially a matrix where you have genes in rows and samples in columns. That allows you to plug in these counts that you get from these tools like RSM or salmon or other quantification methods in order to do the differential expression between the conditions that you have in your experiment. This pipeline doesn't perform any differential expression analysis, and that's intentional because when you start getting involved with stats, that's generally where things start getting complicated. Differential expression, in order to do it properly, you need to factor in all of the various experimental factors you have in your experiment, and there's not really a standardized way of encoding that information.

[7:11](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=431)
To keep things simple, we basically have the RNA-seq pipeline, which gets you to the counts, and then it's up to you how you factor in various sample conditions, like whether you need to factor in the sex of, say, mice, or whether you need to factor in time points in terms of days and how this would affect the differential expression, and other confounding factors that really need to be taken into account. If you want to get an idea of some of the more low-level type mapping types, Reagan gave a great talk last week about the dualrnaseq pipeline, where he explained some of these mapping to various aspects of the genome or the transcriptome and the complications that arise as a result of that. I won't go into much detail there.

[7:54](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=474)
In terms of features, one of the biggest strengths of this pipeline is the fact that it's used so widely. We've got bug fixes. We've got feature requests. We've got contributors coming from all over the world. It's used on various infrastructures and clouds, which, again, is testament to Nextflow itself, and also on various types of input data: small data, large data, medium-sized data, whatever-you-can-imagine type data. That's really one of the biggest strengths of this pipeline.

[8:24](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=504)
In terms of the alignment and quantification routes, we've got three standard ones. We've got STAR and salmon, which Rob Petrow actually helped me add, which was really nice of him to come. He's on nf-core slack. And we went back and forth a bit before I added this functionality. Salmon may not be as widely known, but it also has the ability to take BAM files and quantify from those. And that's the route that we used for the default option in this pipeline. Similarly, there's a STAR and RSM route. RSM is touted to be one of the most accurate quantification methods. And in recent releases, I've really tried to push making this pipeline as accurate as possible to make it a gold standard best practice type pipeline. We've stripped out some of the stuff like feature counts quantification, which doesn't really look at, have any sort of statistical way of modeling where a read count belongs to, for example. There is no feature counts quantification in this pipeline anymore, which is why, actually, HiSAT, you don't have any downstream quantification at the moment, because there isn't an appropriate way to project the reads or the counts onto a transcriptome somehow and then do the quantification, which tend to be the more accurate methods.

[9:35](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=575)
We also have a pseudo-alignment route. So these routes basically skip the BAM file, essentially. They go from a FASTQ file and have this quasi-mapping approach where you use k-mers to then calculate the counts directly from the transcriptome. You skip the BAM file. I guess one of the downsides of that is that it doesn't allow you to get QC of things like genomic contamination and stuff, which you would need a BAM file for. And that's why the major alignment routes at the top here are probably nicer. But there's nothing to say you can't run this and also run this. It's up to you how you run the pipeline. There's an open request for Callisto as well, if anyone wants to help out with that.

[10:15](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=615)
The pipeline runs from bacterial genomes, to all the way to plant genomes, which have ridiculous amounts of duplication. So again, it supports most genomes. There's an inbuilt strand specificity check, which allows you to double check the strand specificity that you've used. This is quite important in RNA-Seq because if you get that wrong, then your quantification will be completely wrong because you're counting reads mapping to the wrong strand, essentially. There's a warning that's currently generated that tells you whether you've got it right or wrong. And a whole bunch of other features like UMI support, RNA removal, genomic contaminant removal, I did recently. And also you can chain this to the nf-core FetchNGS pipeline, which is another pipeline that I've written that allows you to download data just from a set of IDs, SRA IDs, and it generates a sample sheet that you can directly plug into this pipeline. So yeah, various cool features.

[11:12](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=672)
The sample sheet is quite simple. You've got sample FASTQ1, FASTQ2, and strandedness. If you have single end data, you just literally leave out or leave this second column blank and that's it. You have strandedness which I mentioned is quite important for the quantification. There's nothing complicated there. In terms of reference genome options, you only need a FASTQ1 and a GTF or a GFF. If you provide a GFF, this is converted to GTF for the downstream steps. But if you don't provide any, you can also provide indices and stuff to save you having to create them whilst you're running the pipeline. If you don't, then these are automatically created throughout the course of the pipeline.

[11:52](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=712)
There's various parameter docs as well. All of these links work, by the way. I'll make these slides available so you can use them as you go. Genomes, we're looking to move to RefGenie, but the genomes at the moment, we're using Illumina AWS iGenomes. The standard organization is really nice, but it's becoming quite outdated. So we'll be shifting to RefGenie hopefully soon. The results for full-size tests are available on the website. What's awesome about this is that you literally can run a proper full-size experiment with just two parameters. You just need to provide a sample sheet with your samples and the genome and the pipeline will literally generate all of the downstream steps for you. This is available on the website for you to browse. I won't go into much detail here.

[12:37](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=757)
Similarly, there's a bunch of output docs, quite extensive docs about the outputs of the pipeline and some really nice QC plots and stuff that you can have a look at. We're always looking for feedback if we need to improve that. The implementation is nextflow native. It's all DSL2, one process. For each process that we have, we have one biocontainer, and this really is quite modular and it allows us to update and maintain the pipeline a lot easier because each process is essentially its own dependency.

[13:07](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=787)
Nf-core modules, 38 of the modules in this pipeline out of 55 are an nf-core module. Again, it allows us to contribute back to this nf-core modules repository we've created, which is a central repository to host essentially nextflow wrapper scripts for any nf-core pipeline. There's a massive toolkit and stuff that we've built around this to help with maintaining modules and adding them to pipelines. In terms of configuration, one of the most commonly asked questions now with this new syntax is how do I change the process requirements? I've just put some examples down here, but the first thing you would need to do is look in your modules config for the process you want to change. Use exactly the name that is specified in this modules config because it's quite important that you use that because you can have multiple processes with the same name used in the same pipeline if you're using subworkflows and stuff. The logic to select exactly the right process will be already defined in this modules config. Find the process name you want to use. In this case, it's just this that I've copied and pasted out here. And then you can append the arguments as you want. As long as they're non-mandatory, you can append. So here, I just want to add this quality 20 argument. I've created a small config file with these options that will only change the options for this particular process. Similarly, I can change resource requirements if I want, or I can change a container, which is less likely because you want to use a container to ship the pipeline. But if you do, then that's also possible there.

[14:37](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=877)
Differential analysis, as I mentioned before, you get all sorts of counts out that you can use for downstream analysis. The pipeline doesn't do any serious differential analysis. It just generates some basic QC plots for PCAs and heat maps that you can use to straight away figure out how your samples look. But it doesn't actually factor in any sample or experimental information, which you need to take care of downstream. And we're looking for someone that can give us this sort of talk, because it's one of the most commonly asked questions on nf-core, actually, as to what you're doing with the downstream results of this pipeline. And it'd be an awesome bytesize talk to give, actually.

[15:14](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=914)
I've also added this pipeline to Nextflow Tower. So Seqera Labs, which is the home of Nextflow now, and also this product called Nextflow Tower, which is just an awesome way of monitoring and maintaining and administering your Nextflow pipeline executions. We're working hard with the nf-core community, as well, to try and make this even better. There's a community showcase area, the links here, that you can join and get 100 free hours of credits to run this pipeline on Nextflow Tower, amongst others, as well, to show you or to give you a flavor as to what we're doing there. If you want to come and chat with us, you can find us on Slack, and create issues or pull requests on GitHub, on Twitter, and all of these videos and other content is available on YouTube, as well. So thank you for your time. And thank everyone in the Nextflow community and the nf-core community, Biocontainers, and Biocon, and the great infrastructure that they've allowed us to use without reinventing the wheel, and also my awesome colleagues at Seqera Labs. And also, I guess, some of the main contributors to this pipeline in particular, as well, like Mahesh, Gregor, Jose, Phil, who first started it off, and Alex in between, and everyone else that has contributed over time.

[16:36](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=996)
We have a hackathon coming up. If you don't know already, here's a sign-up link I'll put in the slides, but you can find it on the website, as well. The major theme is documentation. If you think we're missing anything, please come and tell us, and we will try and improve documentation wherever we can. Thank you for your time.

[16:54](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=1014)
(host) Thank you, Harshil, for that comprehensive review of the RNA-Seq pipeline. Feel free to ask questions if you have any.

(question) Philip actually had a question. He just wanted clarification on whether we are aligning to the genome in this pipeline, not to the transcriptome.

(answer) Good question. It depends on how you want to look at it. We do align to the genome. You're right. But we project those reads onto the transcriptome, for example, with Rsem, what you end up doing is you get this Rsem, you get this transcriptome BAM, as well as a genome BAM. And the genome BAM is generally what you use for the QC, and the transcriptome BAM is then what Rsem uses to generate the counts. Strictly speaking, yes, we're probably aligning to the genome and then somehow filtering down to then use the transcriptome. That was quite an odd slide. I hope no one noticed, but yeah, we'll utilize them.

[18:04](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=1084)
(question) Okay, and in order of priority of questions, could you clarify how references are built? Do you need a FASTA file?

(answer) Yes. you would need your genome FASTA. You would need some sort of annotation. This pipeline doesn't do any de-novo guided stuff or it doesn't map. It doesn't use just the transcriptome FASTA as an input. If you have a novel species that you've just done, place this on and you've got a transcriptome, but you don't have a proper annotation, this pipeline won't work yet. There's an open feature for that. What the pipeline essentially does is you've got your genome FASTA, you've got your GTF or your annotation, and you extract the transcriptome from those two and use that for all of the downstream analysis. Any indices and any other information is then built from just the FASTA and the GTF.

[19:02](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=1142)
(question) Another question from Philip, he asks, who is Rob Petro?

(answer) Rob Petro is the main author and developer of salmon and a bunch of really other cool tools that are used not only in bulk RNA-Seq, but also now in single cell RNA-Seq for analysis.

[19:26](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=1166)
(question) We have another question from Michael who asks, is it worth considering an R environment with a pre-built DDS object containing all the samples run?

(answer) Sorry, I didn't understand that. Does the pipeline generate a DDS object?

(question cont.) Yeah, like, maybe it will be worth considering an R environment within the pipeline with a pre-built DDS object containing all the samples.

(answer cont.) There is a DDS object, I believe, that is generated at the end of the pipeline for the counts and all of that sort of information. It's a way that you can easily load stuff into your own R environment, but things start getting tricky and then start verging on actually having downstream type analysis like Jupyter Notebooks and all sorts of other RStudio type stuff where you then need to take the results of this pipeline and load them into a more interactive environment. It's something that we've been talking about for quite a while actually, but it's not a trivial thing to figure out, especially when you want to factor in reproducibility and other things and how to do that in a standardized way. It's an interesting question, so at the moment we don't have anything that does it explicitly, but we do generate the DDS file that you can load into your own R environment and do whatever you want with that in terms of the downstream analysis.

[20:46](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=1246)
(question) Okay, so I think you alluded to this before, but Ramon asks, how difficult is it to do conversion from DSL1 to DSL2?

(answer) So for this pipeline it was actually very tricky because it was the first adoption of DSL2 on nf-core and so everything was starting from scratch and I had to basically change things about a gazillion times to actually get to where I wanted to in terms of functionality testing and so on, but now with the awesome infrastructure we've built as a result of various people's learnings over the past year or so, we can now really easily install modules, we've got some really good examples of how to write DSL2 pipelines. I gave a talk about that recently as well, how easy that would be and how you should attempt to tackle it. We can link to that if you follow up on the bytesize channel, I can send you a link to that. It all depends I guess on the complexity of your pipeline, but in theory it should be a lot easier for you than it was for me a year ago.

[21:54](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=1314)
(question) Thank you for that answer and Oliver has a bit of a comment and a question. So he says, great talk, the QC metrics are awesome. Something that will be super helpful for plotting and interrogating QC metrics will be to add the QC results as columns for each row sample, in this case in the input samplesheet.csv. So could the columns in general statistics multi-QC reportable be added to the samplesheet.csv? And he gave an example of how it's done with Tidyverse.

(answer) I mean if you have suggestions as to how we can improve it, we have something similar actually that we've recently just added for viralrecon that I released last week and that's used for SARS-CoV-2 genomic surveillance type stuff where this sort of QC and variant information is quite important. But if you have something functioning already that's even better. If you have an idea as to what we can extract and how we can extract it, if you dump it in an issue and then we can have a look at it. Pull request contributions are always welcome as well. Any suggestions or contributions like that would be more than welcome and I don't see why we can't dump a generic sort of QC flat file type thing, but I think you can export some of that from Multi-QC already.

(answer cont.) Yeah, maybe I can chime in there. So Multi-QC by default will export all tables and quite a lot more into flat files specifically for this reason for downstream analysis. So if you look in your Multi-QC folder, there's the HTML report, but you'll also find the folder called Multi-QC data and inside there, there'll be a whole bunch of files and you can choose what format to have those in as well.

(answer cont.) In fact, that's what I'm parsing for viralrecon. Multi-QC dumps all of these files and it's just really easy not to have to write another parser for every tool that has a log file because Multi-QC is awesome and it does it for you. So I just literally get all of the information from those tables that Multi-QC generates, parse it, and then use that to generate the QC metrics that are reported for viralrecon, for example. I've been meaning to do something similar for RNA-Seq, but I just haven't had the time.

(host) Okay. Thank you, Phil, for chipping in.

[24:13](https://www.youtube.com/watch?v=qMuUt8oVhHw&t=1453)
(host) I don't know if there's anyone who has a question, would like to unmute, but if it isn't the case... Thank you, Harshil, for this splendid review and answering the questions quite well. I guess we'll see each other. We'll see everyone in the next bytesize talk next Tuesday.

</details>
