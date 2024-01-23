---
title: 'Bytesize: nf-core/nascent'
subtitle: Edmund Miller - University of Texas at Dallas, USA
type: talk
startDate: '2022-11-01'
startTime: '13:00+01:00'
endDate: '2022-11-01'
endTime: '13:30+01:00'
youtube_embed: https://www.youtube.com/watch?v=chayGGPTnfM
locationURL:
  - https://doi.org/10.6084/m9.figshare.21444867.v1
  - https://www.youtube.com/watch?v=chayGGPTnfM
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: nf-core/nascent

This week, Edmund Miller ([@Emiller88](https://github.com/Emiller88)) will talk about the newest developments in the nf-core/nascent pipeline.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=chayGGPTnfM&t=1)
(host) Hello, everyone. My name is Franziska Bonath. I'm the host for today. With us is Edmund Miller. He is a PhD student at the University of Texas at Dallas. Edmund is going to talk today about the pipeline from nf-core called nascent. And to you.

[0:19](https://www.youtube.com/watch?v=chayGGPTnfM&t=19)
Hey, good morning, everyone. I'm Edmund. Let's get started. If I can advance the slides. Okay, a quick overview of what we're going to talk about today. A quick background on nascent transcript identification. Because I'm not sure if it's as common as some other assays such as ChIP-seq, RNA-seq. A brief history of the development of the pipeline. Lastly, we're going to talk about the pipeline itself and give a brief overview of that.

[0:55](https://www.youtube.com/watch?v=chayGGPTnfM&t=55)
Quick background on nascent transcript identification. The goal is to identify the changes in transcription of the RNA and what's going on in the cell at that specific point in time. Rather than, say, RNA-seq, which isolates all of the RNA in the cell at a steady state. That would be like your mRNA and things that have matured versus what's actually being transcribed. You can get an actual response to things like heat shock or viral infection. Pulling out the transcription activity sites through metabolic labeling of these. We won't go into that too much today, but I'm happy to discuss that with anyone in the future. The problem with that is that we're covering a lot of different assays, not just one or maybe some slight variations with that. We're covering lots of different chemistries, lots of different steps, et cetera. Because of that, some slight variation in the computational pipeline can actually lead to 25% change in the results of the transcript calling. We'll talk about why that is later.

[2:05](https://www.youtube.com/watch?v=chayGGPTnfM&t=125)
Specifically, what I'm interested in are enhancers. There's a lot of different things that you can call with nascent transcripts, such as miRNAs and long non-coding RNAs. You can call the gene sequences as well, but specifically, what you can pickup that isn't really possible with RNA-Seq, is enhancers. These are cis-acting DNA sequences that can then increase the transcription of genes. A lot of people are probably familiar with promotors, they (enhancers) act in tandem with promoters. Part of the problem with these enhancers and identifying them is that there's hundreds of thousands of enhancers, but we have this evidence that the enhancer promoters interact through various other assays such as 3C. We also have evidence from these nascent transcript assays that enhancer RNAs are produced at these enhancers. They have a very short half-life, and they're in low abundance. We don't usually pick them up in general RNA-Seq.

[3:10](https://www.youtube.com/watch?v=chayGGPTnfM&t=190)
This is just a quick infographic of what's happening here in the enhancer promoter looping. Over here on the right, you can see that we have the promoter and Pol II. Then, we have the mRNA coming off. This is what everyone's probably very familiar with. It's being produced. This is what you pick up in RNA-Seq. Then, you also have transcription factors and cofactors, but what we're really interested in, or I am specifically, is this enhancer on the other side and the eRNAs coming off of that with the Pol II activity. Those are thought to pull in all these transcription factors and cofactors and all these various other things.

[3:50](https://www.youtube.com/watch?v=chayGGPTnfM&t=230)
What do the reads actually look like since we're talking about bioinformatics here and what we're interested in? These are just a couple of the various assays. This is out of a recent paper that I thought was really good that summarizes all these. You can see we have GRO-cap, csRNA-Seq, NET-CAGE, STRIPE-Seq, PRO-Seq, all these. But let's start here at the bottom with total RNA-Seq. As you can see, just to orient everyone, we have the known enhancer here in yellow. Then, we have the reads over here on the left along the gene that's known. In total RNA-Seq, there is a strong peak on the antisense. You can see that we go along and have some reads coming from there.

[4:40](https://www.youtube.com/watch?v=chayGGPTnfM&t=280)
The main point here is that we don't pick up the known enhancer in total RNA-Seq. There's just not enough reads and not enough mature RNAs happening. Whereas something in GRO-cap, for example, you can see that we really pick up the known enhancer and have a lot of signal coming from there. However, we don't pick up the entirety of the transcript in GRO-cap, for example. But you can see we also pick up this opposite transcriptional start site that's going in the other direction from the gene body itself.

[5:16](https://www.youtube.com/watch?v=chayGGPTnfM&t=316)
There's other things like PRO-seq, which actually are nascent transcript assays, where you can see we pick up a little bit of the known enhancer. We don't have such a pronounced peak, perhaps. But then, we also pick up the entirety of the gene body and things that are being transcribed all the way along.

[5:31](https://www.youtube.com/watch?v=chayGGPTnfM&t=331)
As I was just talking about, we have two different kinds of assays that we're supporting. We have nascent transcripts, and then we have transcriptional start sites. I think this image from the same paper does a great job of illustrating this as well. Part of the problem is there's like 13+ assays for nascent transcript identification and transcriptional start sites. And as I said before, minor changes in the sample processing could lead up to greater than 20% in the final results. That's what they found, and I was validated by them. I'll talk a little bit more about that in the history. So let's start down here at the bottom. You can see the promoter and the true transcriptional start site here. The blue is the TSS assay, like GRO-cap, that I was just talking about, whereas the nascent transcript assay would be PRO-seq. Neither of these are generic RNAs either of what you're thinking about. So you can see in the TSS assay, we get a very pronounced peak and this is actually at the promoter sequence at the very beginning of the promoter. Then we have a slight break, and this is CpG island here. Then you have the nascent transcript assay. That picks up the entirety of the gene body and the elongation of that. So these are the two different types of assays that we're picking up.

[6:53](https://www.youtube.com/watch?v=chayGGPTnfM&t=413)
The interesting part is that we're picking up enhancers as well, based on those. We have a TSS assay, and that's where we're picking up the initial transcription start site. Then we can also pick up the entirety of this transcript and where Pol II is actually working along the entirety of it. Over here, this is just talking about the directionality of these and whether we're pulling them with a cap or not a cap. I highly recommend the paper if you're interested in that.

[7:24](https://www.youtube.com/watch?v=chayGGPTnfM&t=444)
A quick history of the development: version 1.0 was developed by Ignacia Tripodi and Margaret Gruca, and was released April 16th in 2019. In Parallel in 2017, the Tae Hoon Kim lab at UTD started working to reproduce a paper that came out in 2018 in a second data set, and I was mostly responsible for that. This is where I got my start with bioinformatics and reproducible research because I struggled to build a reproducible pipeline and reproduce the results from that paper. That's where I kept getting into the 20% variance of these things can really make or break the transcript calling. I didn't understand that at the time, but now after being validated, it feels great that it's so much different than some other assays where the bioinformatics pipeline doesn't affect it that much.

[8:23](https://www.youtube.com/watch?v=chayGGPTnfM&t=503)
I started creating my own CI/CD workflows and templates for SnakeMake in around January 2020. As soon as we had a little lab hackathon, introducing it to everybody, I found nf-core the week before and started looking to move everything over to that because I was excited to work with others on that and doing a lot of great work here.

[8:44](https://www.youtube.com/watch?v=chayGGPTnfM&t=524)
So let's talk a little bit about the pipeline. This is how far we've come. This was a SnakeMake DAG because there wasn't a DAG of the V1 of the nf-core pipeline, but this is what I had in 2018. You can see the original presentation where I'm boring my lab with things like Docker and other things as well in that, but you can see the majority of this is we're just handling Homer and alignment, it's pretty much all, and then maybe an intersection of histones and the GM data as well and handling those two cell lines, very rudimentary.

[9:25](https://www.youtube.com/watch?v=chayGGPTnfM&t=565)
This is the obligatory metro map that I finished last night, and then James already has some feedback for me, but I like to thank all of those who worked on that, it was a great template and really easy to get going with that. Let's start over here with the .fastq, and then we can pretty much zoom through everything here because we're really standing on the shoulder of giants here and using a lot from RNA-Seq, which is great because it's a much smaller use of pipeline and there's a lot less users, but we benefit from all of those bug reports now, with subworkflows and modules and all those other things. We can really jump all the way to transcript identification. We just make some genome maps up here, these are the only unique thing to us from RNA-Seq, and we support a few different aligners.

[10:18](https://www.youtube.com/watch?v=chayGGPTnfM&t=618)
The first thing is we're grouping all the replicates up. Basically that's anything that's a technical replicate that we want to group up, to increase the signal and biological replicates. Then we feed that into for GRO-seq over here, if they're specifically, because that's what I've been so interested in. We feed that into Homer and GroHMM optionally, and we'll talk a little bit more about that. Everything else that's a transcriptional start site and GRO-seq and others. We feed that into PINTS as well, and then we go into bed tools and we can intersect the two of these with a filter and without a filter, and then basically only call regions that we're interested in and drop the regions that we're not. We can drop the regions that are gene bodies and promoters because we know that those aren't going to be eRNAs or other interesting RNAs. We can also make sure that we keep only regions that we're interested in, such as like those with histone modifications that indicate eRNAs. Then we just do some quick quantification, and then we move into MultiQC.

[11:29](https://www.youtube.com/watch?v=chayGGPTnfM&t=689)
Another little added benefit that we were interested in was supporting CHM13, which is a new reference genome that came out recently. Highly recommend y'all look into that as well, if you're interested in that. I'll be adding this to the template soon. But the main thing here in this infographic that they found is they were specifically looking at methylation data and how the new reference improved calls.

[11:56](https://www.youtube.com/watch?v=chayGGPTnfM&t=716)
You can see over here on the left is the number of MACS peaks, and then you can see the blue is the old reference. The CHM13 reference are the additional calls that were made from using this reference. So these may not be much and may not be of interest in things that are well known and well understood, but very relevant for nascent transcript calling. We have support for that in our IGMs config, and you can just use CHM13 and align to that.

[12:29](https://www.youtube.com/watch?v=chayGGPTnfM&t=749)
Let's talk a little bit about the transcript identification because that's the most interesting part of the pipeline and what makes it unique. There's a couple of different options, as I said. First, if you're doing GRO-seq, I have some great support for that. If anyone would like to support other assays or would like to see it supported, please open an issue. So first is GroHMM, and this is what kind of sparked us getting into the nextflow and going into bigger pipelines. It was difficult reproducing this and running this on big enough machines to actually use it because it's an R package. It was released in 2015 by Minho Chae, Charles Danko, and Lee Kraus, actually just down the street at UT Southwestern.

[13:20](https://www.youtube.com/watch?v=chayGGPTnfM&t=800)
As you can see by the graphic up here, GroHMM greatly outperforms HOMER in just about all of these metrics. SICER is actually just a chip-seq calling or an old chip-seq peak calling algorithm. It actually outperformed HOMER, which we thought was interesting looking at this graphic. There's a couple of drawbacks to GroHMM though. It's very time consuming because it requires tuning and it's also quite memory hungry when you're running on a bunch of samples. We also reached out to the authors and Charles Danko recommended that we use T-units, which is an unpublished R package that doesn't require tuning. So stay tuned on that. But right now, GroHMM works and it does perform very well.

[14:06](https://www.youtube.com/watch?v=chayGGPTnfM&t=846)
This is calling the entirety of the transcript though, just to note up there on the left... oh... I think I missed HOMER. So I'll just talk a little bit about HOMER then. Without it (the slides). It uses a little bit more naive of a peak calling method. It's just looking for the transcript and the difference in the peak and itself. On those, it was released in 2010 out of the glass lab, and it was maintained by Chris Brenner. That was what we originally used in our paper. It works pretty well. The problem is it was made in a land before Docker. So it has a couple of problems with the way that it wants to pull in the references for you. But I finally realized you can just pass a FASTA in. It's like one line in the documentation and that works amazing. We just run that on everything because if you're going to run GroHMM and wait a couple hours, 20 minutes with HOMER, you might as well get some results on that as well. So again, I missed the slide on that.

[15:18](https://www.youtube.com/watch?v=chayGGPTnfM&t=918)
Let's now jump into PINTS identification. This is a new assay that just came out in 2022, and it's very exciting. I just left this in up at the top in Figure A. This is just also illustrating the difficulty in reproducing these, and this is on the exact same data sets. You can see the difference in the HOMER and the GroHMM results and just how much they vary by just a slight tweaking of a tool. You'd expect maybe better performance, but you wouldn't expect completely different results based on what tool you're using. Down here at the bottom, this is the PINTS identification method. Just in a rudimentary way, it works very similarly to MACS2. This is a potential true peak based on the density of this, and it's very easy to pick out. It does some algorithms, picks up the local background noise from these, and these are the light blue. You can see in the purple from those, that's then a potential peak that it needs to test and see, is that actually a peak or is it just more noise from the assay itself? What PINTS is doing is really picking up these TSS start sites. As you can see from the read pile up here, it's just picking up the TSS site rather than the entirety of the transcript, which might actually lead out all the way along here.

[16:45](https://www.youtube.com/watch?v=chayGGPTnfM&t=1005)
As I said, it was released in 2022. So it's a little more relevant than 2010 and 2015. from the the Yu and Lis lab. It determines the TSS start site is really what it's doing as opposed to the entire transcriptional unit because it's mainly focusing on TSS assays. It also achieves the optimum balance among - this is from their paper - resolution, robustness, sensitivity, specificity, and computational resources required.

[17:17](https://www.youtube.com/watch?v=chayGGPTnfM&t=1057)
There's a couple of other tools that can also be used such as D-reg, but those required GPUs, and you start getting into all kinds of difficulty for users and specific machinery. It also supports TIN assays just out of the box and works. So that's a quick win, and then we can kind of support all of those through using PINTS and just handling most of the upstream and downstream processing of those. Cunningham's Law here, the best way to get the right answer on the internet is not to ask a question, it's to post the wrong answer. If you think that any of this information isn't correct, or we should be doing things differently, please open an issue or drop into Slack. I know there's not a lot of cohesion on the nascent assay transcript identification, but I'd love to help the community build a kind of a group, ideal workflow on this. So with that, I'll take any questions.

[18:24](https://www.youtube.com/watch?v=chayGGPTnfM&t=1104)
(host) Thank you very much. I have now allowed everyone to unmute themselves. If there are any questions, you can... yes, Harshil.
(question) Hi Edmund, yeah, thanks, great talk. I think, for this pipeline in particular, as we're realizing now a lot on nf-core is that we've got the really nice pipelines, but we need to be able to validate the results between releases and stuff. This is the thing that have come up during the summit. I think this is a really nice example of that, because as you mentioned, you tweak some parameters, or you run the pipeline in a different way, and you get all of that variability in the results. It's really important to be able to reproduce the results. Have you thought about full size test data sets and how we can validate whether the results are actually optimal across releases? So say, for example, you or someone else comes to tweak the pipeline, that we're not negatively impacting the results that you should be getting.
(answer) Exactly, that is something that I've thought about. I haven't gotten a AWS full test going for GRO-seq yet. I do have two tests that were in the PINTS that they created some test data examples that I asked for, because they didn't have any examples of the actual usage of it. From those, we can then call the peaks on CoPro and the other ones. I'm missing the other one. But there's two test data sets already that are full data sets that I ran. Then I have regression tests of those that I'm saving as well to compare against. They actually have an entire element matrix. We can probably pick a few of those and see if we can reproduce those each time, or at least benchmark where the nascent pipeline is and make sure that we're not changing drastically on those.

[20:30](https://www.youtube.com/watch?v=chayGGPTnfM&t=1230)
(question) That would be awesome. A second question. So no controls, right? You don't have controls for GRO-seq?
(answer) The control sample is kind of included into that. They talk about in the PINTS paper a little, some tools require you to have controls. The tools that we're using don't really require controls.

(question continued) Okay, so then the background model is built up and then the caller will call the peaks based on some random distribution in the genome.
(answer) Yep.

[21:00](https://www.youtube.com/watch?v=chayGGPTnfM&t=1260)
(question) Last question, why not using MACS and other conventional callers? Why use Homer? Homer seems quite primitive, I guess, in terms of peak calling and stuff. Why not something more sophisticated like MACS? Is there more false positives?
(answer) Legacy. Homer, you can also tweak some of the important things of like it picks up on the... I missed the image, but basically it picks up on the peak and then it picks up on the trailing tail of it. That is actually the piece that's really important there instead of... Here, I'll just pull it up. This is what Homer's actually doing. Whereas in MACS, you might just pick up the peak. You're actually picking up this downstream transcript is why Homer's unique to that.

[21:46](https://www.youtube.com/watch?v=chayGGPTnfM&t=1306)
(question) Okay. SIZER presumably does something similar because it calls larger peaks as well, right? It's able to call these sorts of counts?
(answer) uhum.
(question continued) Okay, cool. Thanks a lot, man.

[21:57](https://www.youtube.com/watch?v=chayGGPTnfM&t=1317)
(host, question) There's also another question in the chat. Why do you use feature counts and not other quantification methods as in RNA-Seq?
(answer) Feature counts is always just what I've used for that. I'm open to other ideas on it. It's not the exact same as RNA-Seq and most of those are RNA-Seq specific, is part of the issue on the quantification of those. So the difference is we pass in the genes, count with those. Then we also count with the identified transcripts and identified transcriptional start sites of those and give you counts of all of those. That's the difference. Downstream you have to do your own math behind the scenes and stats because it's not the exact same as RNA-Seq in terms of how the math works out on those. Again, also not well-defined.
(audience) You're counting with RNA-Seq, you're counting things that overlap, spliced transcripts, whoever's GRO-seq, you're looking at the entire gene body where splicing isn't important. Feature counts can do that in this case, whereas with RNA-Seq, as we've known and had previous discussions, it's not ideal for the transcript splicing type quantification.
(answer continued) Exactly. Exactly. Well said. It's just... it can work in a very simple way is the reason that we're using feature counts.

[23:27](https://www.youtube.com/watch?v=chayGGPTnfM&t=1407)
(host) Okay. Thank you. I don't see any more questions. So with that, I want to thank you, of course, Edmund, but also the Chan Zuckerberg Initiative for funding the bytesize talks and as usual, if there are any questions, you can always go to the nf-core workspace on Slack and the nascent channel and ask your questions there. Thank you very much.

</details>
