---
title: 'Bytesize: nf-core/rnafusion'
subtitle: Annick Renevey - Karolinska Institutet
type: talk
start_date: '2022-09-13'
start_time: '13:00 CEST'
end_date: '2022-09-13'
end_time: '13:30 CEST'
embed_at: 'rnafusion'
youtube_embed: https://www.youtube.com/watch?v=iP47pokiPB4
location_url:
  - https://www.youtube.com/watch?v=iP47pokiPB4
  - https://doi.org/10.6084/m9.figshare.21206537.v1
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: nf-core/rnafusion

This week, Annick Renevey ([@rannick](https://github.com/rannick)) will introduce the nf-core/rnafusion pipeline!

<details markdown="1"><summary>Video transcription</summary>
**Note: The content has been edited for reader-friendliness**

[0:01](https://www.youtube.com/watch?v=iP47pokiPB4&t=1)
(Maxime) Hello everyone, Maxim here. I'd like to welcome Annick Renevey from Clinical Genomics, Karolinska Instituted. She's going to talk about the rnausion pipeline, which is a pipeline I personally like a lot, because I used to work on this one. I'm still helping out a bit, but she's doing a way better job tham I was doing at the time. Thank you Annick for that already!

[0:30](https://www.youtube.com/watch?v=iP47pokiPB4&t=30)
(Maxime) Before we start I just like to quickly thank the Chan Zuckerberg Initiative, for helping us in organizing this bytesize talks and to all of our listeners. You'll be able to unmute yourself at the end of the talk for questions.

[0:47](https://www.youtube.com/watch?v=iP47pokiPB4&t=47)
Hi everyone! I'm Annick, I'm the main developer currently of rnafusion and I will lead you to an hopefully short introduction what is our goal with the pipeline.What we can get out and and how to use it in a few a few words. I will start with our angle of the pipeline. We are coming from the Clinical Genomics unit in Stockholm. It has a lot of links with clinical diagnostic and we are providing analysis tools that help diagnosticians when they are reporting back to patients. So we are really part of rooting clinical care. Clinical Genomics is sitting in SciLifeLab, which is a conglomerate of four different universities, that work together, so we are part of an organization that is part of an organisation, hence the multiple applications.

[2:01](https://www.youtube.com/watch?v=iP47pokiPB4&t=121)
We need fusions because they have been detected increasingly in in many common cancer types and they are a very valuable tool for Diagnostic purposes. The first versions were developed during Martin Proks' Master Thesis at SciLifeLab and it he has been building on the work of others. Maxime has been contributing to it a lot, but we have a lot of of different contributors along the years. Unfortunately it got outdated - as many scientific software do - because we got a lot of other things on our desk. The software got updated, the database got much better and and all of a sudden, when we wanted to use it end of last year, the pipeline was effectively broken. We couldn't download the references that we needed to run, so there was need for some rework.

[3:07](https://www.youtube.com/watch?v=iP47pokiPB4&t=187)
That's when I came into play. There is the now version 2.0.0 that has gone out already, it was a complete rewrite and an upgrade to the DSL2 syntax. It includes flexibility, so that you can make the pipeline do more or less what you want. Otherwise just open an issue and we'll see what you can do with CLI options and adding visualization and quality control tools.

[3:46](https://www.youtube.com/watch?v=iP47pokiPB4&t=226)
The main goals is to detect Fusion in RNA sequencing, but there are many different ways, different tools, to detect fusions. The idea is to combine the power of the tools available and to compare them. To compare them between themselves and also with databases or fusions that are already present This can help you in case you're looking for a common Fusion type but. if you're looking for a novel Fusion you might want to go further than just a database, so this is just an indication. The pipeline is also completed with visualization tools and quality control, so the pipeline overview looks like this.
You can imagine it like a network of different Subway Lines. You can take any of the Subway Lines, all of them or just maybe Ariba, SQUID and pizzly and maybe you don't care about FusionCatcher and STAR-Fusion. You also have a parallel line that is consisting of the quality control and the core analysis tool, which will lie here, where the fusion report (which is a tool developed by Martin Proks) that basically takes all of the fusions detected by the five different software lines, put them together and checks if this Fusion is identified by this tool and is present in this database. Once we have looked at this, we take every Fusion that has been identified by two tools or more and we look again in more in detail into it with Fusion inspector, collect that statistics Etc

[5:42](https://www.youtube.com/watch?v=iP47pokiPB4&t=342)
Here is how the output of fusion report looks like. You can see that it has an interesting dashboard where you have all the tools, known versus unknown, Fusion databases and by how many tools the Fusion was detected. In our case the tool was very sensitive, so it detected many Fusion. Now I'm gonna try to do an interactive demo. Let's see how it works. The table here is very nice to look at it a bit more in detail. You can see that here I can highlight how many fusions were identified with pizzly and if I hide the fusions identified by pizzly, I can have a look which tool identified how many fusions. This is quite interesting. If I remove fusions that were detected by one tool probably pizzly, then I have a bit more of a detailed panel. This table though is is very interesting because you can sort how you want it, like change the orders Etc.

[7:30](https://www.youtube.com/watch?v=iP47pokiPB4&t=450)
This is a sample that is artificial, you will probably never see this - hopefully - in an in a natural sample. This is a sample consisting of 20 fusions so you have those in the sample and as you can see, they are found out by all the tools corresponding to the 5 tool hits. You have a scoring function that depends on the number of tools that have found the fusion and also on the different databases that have found it.

[8:14](https://www.youtube.com/watch?v=iP47pokiPB4&t=494)
That's that's a quite valuable tool if you want to compare between different tools. Now coming back to the different results we can have a look at Fusion inspector. Here is just one side the HTML output of Fusion inspector. There is a lot more and I really encourage you to run it and look for yourself. There is something that can be of interest to you. If you're interested in a special part of the fusion there are bam files, there are a lot of tables of Statistics so this is just an overview.

[9:02](https://www.youtube.com/watch?v=iP47pokiPB4&t=542)
It's an interactive table, so you can again look at at fusions and you have also some visualization possible in the browser, among others. You can see some statistics here and a bit more about the Gene and and their positions.

[9:30](https://www.youtube.com/watch?v=iP47pokiPB4&t=570)
The last visualization tool that I wanted to show you now is the Arriba visualization tool. It's only done for fusions that have been identified with Arriba. You get a PDF file out: one slide per fusion. This is one Fusion. You can see a very detailed view of the breakpoint. You can really have an idea of the sequence Etc. You can also have a quick look at the retained protein domains which might be important pathologically and a few supporting read counts, like statistics.

[10:26](https://www.youtube.com/watch?v=iP47pokiPB4&t=626)
About how to use the pipeline. What you would have to do is first build the references. This requires patience because at the time we are building the STAR-Fusion reference from scratch and that takes about 24 hours on an HPC. Don't be surprised that it takes a long time, it is what it is for the moment. I'm hoping to make it shorter at some point. If I can host the build references directly, but this is something I'm working on. You would have to start by creating a COSMIC account and passing your your username and password to the software. Then you would have to specify --build-reference references and the tools that you want. I put "all", because I find that it makes sense to build for all tools, but if you only want to use Arriba then you can just do --arriba. If you want to use Arriba and fusion capture you would choose to --arriba --fusioncapture, and you would only download Source references.

[11:58](https://www.youtube.com/watch?v=iP47pokiPB4&t=718)
Then you need to provide genome space which is a pass to your references, and --outdir, which will be the output directory of the run. In this case it will not contain very much, because all of the data, the references, will be generated in genome space. You will have still the execution Tracer logs, the versions Etc in the outdir. If you don't specify --build_references it will run the actual analysis. You have the possibility to do all of the analysis, or you can just do any combination of the four tools that you want. If you want Fusion capture on SQUID --fusioncatcher, --squid. If you want everything except pizzly, just specify each tool and not pizzly.

[13:05](https://www.youtube.com/watch?v=iP47pokiPB4&t=785)
You also need an input this time. It will not complain if you do not have an input in when you build the references, but if you try to run the pipeline it will complain if you don't have a sample sheet. You need to create a sample sheet with your sample. The first three columns are standard nf-core: a sample name, fastq_1, fastq_2 and on top you have the strandedness which depends on your library preparation kit. You need to link the genome space, passed to your references, and the outdir is this time very important, because it will contain all of your analysis.

[13:48](https://www.youtube.com/watch?v=iP47pokiPB4&t=828)
I included a few things to help you gain more flexibility in your usage of the pipeline. You might just use it very standard, you don't have to even look into these options. But if you are looking into doing something more specific, or gaining some time at runtime, then it might be useful for you. You have the possibility to skip the visualization, if you're just interested in the different results for the tools but not the visualization. With skip_vis you will skip Arriba visualization and fusion inspector. skip_qc will skip the entire QC line. You could manually feed references paths for each tool if you have them in different directories and you can also just run Fusion inspector with the option fusioninspector_only and then you will have to provide Fusion inspector fusions and the paths to a file that you manually construct, and that has a fusion that you want to to look into for this sample. Then only Fusion inspector would be run.

[15:06](https://www.youtube.com/watch?v=iP47pokiPB4&t=906)
As you can see you can have a few possibility to enter at different points of the pipeline. It has been suggested to me to maybe also add in alignment shortcut, so you would feed manually alignments to the pipeline. This is a great idea. At the moment, as you can see, we are aligning basically for each line. This is because we have each time parameters for the alignments that are optimized for the different Fusion detection tools. It should perform slightly better, but if you want to save time you might want to bypass the steps. This is something that I think I will work with in the very near future.

[15:56](https://www.youtube.com/watch?v=iP47pokiPB4&t=956)
Speaking about the future, I will talk about what's going on. There will be a next release, hopefully very soon, with trimming: adapter and quality trimming. With the possibility to run stringtie as an extra line. That will be helpful because there is a type of fusion that is not detected by any tool currently implemented and that is when you skip for example an exon if it was in the same gene. That should be resolved by using stringtie. But again if you're not interested in this type of fusion, you could skip it completely or run just this.

[16:44](https://www.youtube.com/watch?v=iP47pokiPB4&t=1004)
I'm really looking forward to see if we find a solution for the AWS Mega test, so that we can host a demo on results on the website. Then you can have a look yourself at the results that the pipeline can give you.
On the in Nextflow summit in Barcelona I will present more details, hopefully on our implementation in production: what sort of issues we are facing about data and I'm hoping also to release a how-to video with more details and hands-on demonstration about each command line option.

[17:31](https://www.youtube.com/watch?v=iP47pokiPB4&t=1051)
If you have any questions I am happy to hear them now, or feel free to reach out to us on slack or on GitHub. open an issue. It's great to hear about the different experiences. Thank you for for your attention.
(Maxime) Thanks a lot Annick, that was super good and clear. I really like it. Now people should be able to unmute themselves if they have any questions.

[18:18](https://www.youtube.com/watch?v=iP47pokiPB4&t=1098)
(Question) Thank you, nice talk! The visualization and everything is amazing because it really puts these things into perspective. For a while this pipeline stood out here, because we don't really have that sort of thing for most other nf-core pipelines (other than maybe the multi-qc report). So it's quite nice having something that's customized, that you can use to organize and query the results and stuff. I guess we should probably start thinking about how we do that for other pipelines as well at some point.
In terms of references though, is it just human samples that the pipeline works with or do other sample types work? I tend to keep up to date with what's going on in the slack Channel, but things get out very quickly on nf-core as you know. What I've always been confused about is about the references that are out of the box compatible with the pipeline, and how easy it is to create references to use with the pipeline. That's been quite a big issue recently, in terms of creating these references and using them.
(Answer) I won't say it's easy but it's possible. You could basically feed any reference that you build yourself to the pipeline. The recipes are there in the pipeline when I built them. If you feed an non-human - mouse or something like that - fastq and gtf, then you would be able to build the references for non-human. Not guaranteeing that it's easy.
(Question) so what's the problem, where's the complication?
(Answer) you might not find the exact same types of files or you might be missing databases. The whole database is human based, so you won't be able to compare. Mostly I think it should be possible. It's more about searching for the right files, testing a bit,...
(Question) ok. Do we have these databases on nf-core anywhere? That are just easily pullable? Or is that part of what you were going to do next?
(Answer) That's something that I would love to have and that would reduce our tests a lot. At the moment, as I said, it takes 24 hours to build the references. If I could host them somewhere, if you have some space, shout out!
(Question) I'm sure we could make it available somewhere. I guess, going slightly off topic now, but the AWS iGenomes bucket we have an S3 has typically been used just for that and so we haven't really added any other custom files to it. It could be something we could maybe just push there, if it's tested and it works and we know. Then we can just have a sweet path and upload it. If you get a list of assets together, then I think maybe maybe we can host in iGenomes and we can try and make it happen. Thanks for the talk and see you in at Summit!

[21:54](https://www.youtube.com/watch?v=iP47pokiPB4&t=1314)
Good, anyone has any other question? Then I guess we are good. Thank you very much again Annick.

</details>
