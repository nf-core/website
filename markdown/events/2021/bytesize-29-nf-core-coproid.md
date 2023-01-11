---
title: 'Bytesize 29: nf-core/coproid'
subtitle: Maxime Borry - MPI-SHH, Germany
type: talk
start_date: '2021-11-30'
start_time: '13:00 CET'
end_date: '2021-11-30'
end_time: '13:30 CET'
embed_at: 'coproid'
youtube_embed: https://youtu.be/gU4jx1pb8Tw
location_url:
  - https://youtu.be/gU4jx1pb8Tw
  - https://doi.org/10.6084/m9.figshare.17099534.v1
  - https://www.bilibili.com/video/BV15b4y1B7Uc
---

# nf-core/bytesize

Join us for a special pipeline-focussed episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 29: nf-core/coproid

This week, Maxime Borry ([@maxibor](https://github.com/maxibor)) will tell us all about the nf-core/coproid pipeline.

<details markdown="1"><summary>Video transcription</summary>
**Note: The content has been edited for reader-friendliness**

[0:01](https://youtu.be/gU4jx1pb8Tw&t=1)
Hi, everyone. As usual, I'd like to begin by thanking you for joining us and the Chan Zuckerberg Initiative for funding all nf-core events. We're joined today by Maxime Borry from the Max Planck Institute for the Science of Human History in Germany, and he will be presenting the nf-core/coproid pipeline. CoproID has been described to help identify the true maker of Illumina-sequenced coprolites or paleo feces by checking the microbiome composition and the endogenous DNA. I'm curious to learn how that's done. I am so excited for your talk today, Maxime. If you have any questions for Maxime, you can either unmute yourself at the end of the talk or use the chat function, and I will relay the questions over to him. Thanks very much for agreeing to present for us today, Maxime. I'd like to hand over to you now. Over to you.

[1:06](https://youtu.be/gU4jx1pb8Tw&t=66)
Okay. Thank you very much, Renuka. Hi. As Renuka mentioned, I'm Maxime. I'm a doctoral researcher at the Max Planck Institute for Evolutionary Anthropology. We moved from an institute this summer. We're now based in Leipzig in Germany. Let me tell you in the next 15 minutes about this coproID pipeline that I developed and published last year. What we're going to talk about today. We're going to briefly talk about what is a coprolite, the challenge of identifying their source or sometimes their sources. The solution that we brought with coproID, and also, because coproID was published more than a year ago, I will briefly touch at the end about nf-core then and now.

[2:03](https://youtu.be/gU4jx1pb8Tw&t=123)
What is a coprolite? If you google coprolite, this is probably the picture you're going to end up on. This coprolite is actually quite famous. It even has its own Wikipedia page. It's known as the Lloyd's Bank coprolite because it was discovered when they were doing construction work at the Lloyd's Bank, I believe in London. It's from the 9th century, so a bit more than 1000 years old. This coprolite is now on display in a museum. The fun fact about it is, one day while they were visitors, it broke in three pieces. It was someone's job to re-glue this old poop back into one single piece. If you read in the archeological literature, you will very often find the two words coprolite and paleofeces used interchangeably. They are not exactly the same thing, but for the purpose of this presentation and most of the time they will be used interchangeably.

[3:10](https://youtu.be/gU4jx1pb8Tw&t=190)
Why do we study coprolites? Coprolites are very interesting because they're in the window into the past of the gut microbiome of ancient individuals. It's a way of studying the gut microbiome as they are in modern studies, but for ancient archeological samples. For example, there was this paper that was published I think a month ago where they looked at coprolite from different time periods, including one from the Iron Age. 2,500 years ago, and they found proofs of blue cheese and beer consumption in miners from Austria. Without coprolite, you wouldn't be able to prove that they were consuming this blue cheese and this beer. It's quite cool. You can also look at things such as diseases and a lot of other things.

[4:14](https://youtu.be/gU4jx1pb8Tw&t=254)
When you work with ancient poop samples, there is an additional challenge of identifying the origin. You work with a modern sample, this question is quite straightforward because you sample directly from the origin. Let's say you do a study of the gut microbiome of Elmo that you're going to eventually successfully publish in nature. You know who did the poop that you're going to sequence. You know whose microbiome it is because you're directly getting the sample from the source. But when you work with archeological samples, it's not so straightforward. Coprolites are often found in the archeological context where you can't really easily attribute them a maker. You can't really say who made the coprolite either because there are no nearby human remains. Coprolites are very often found in isolation. For example, at the bottom of a mine, you know that there was human activity, but you don't have any skeletons, so you don't know really who made it for sure. Sometimes you're more lucky. For example, there were coprolites that were found directly in the guts of mummies. In this way, there is no question, but some other times it's not so obvious.

[5:44](https://youtu.be/gU4jx1pb8Tw&t=344)
The shape and the consistency is quite often degraded. For example, below you can see the picture of coprolites that we used in the coprolite article. If the archeologist didn't identify them first as poop, I wouldn't have even guessed that they were ancient poop. Very often in ancient sites, you have on the same place people were living with their animals, especially pet animals that were domesticated animals, like dogs and pigs. You had mixed human and animal occupation, meaning that the author of the poop could be different possibilities.

[6:31](https://youtu.be/gU4jx1pb8Tw&t=391)
What we came up with in the coprolite identification pipeline is the following. After some pre-processing, we go into two different parts, two different ways of identifying the host or the maker of the coprolite or the paleofeces. The first way is by comparing to the reference genomes, genomes with an "s", I'll come back to it in a second. The second way is to do some metagenomic profiling with some machine learning to identify the host from the microbiome composition. Then at the end, we put them together to give a nice report to the user.

[7:15](https://youtu.be/gU4jx1pb8Tw&t=435)
When you do host DNA competitive mapping, you align the DNA sequence that you have against your most likely genomes. For example, here, I took a modern pig microbiome study, and I aligned the reads that they got to the human and the pig genome. I looked at the log-fold to change of human versus pig or the pig versus human. Because it's a pig microbiome study, you're expecting to find much more pig DNA than human DNA in your sample. For most of the samples, it's true. There is, however, one surprising sample, and I let you come up with an explanation by yourself, but normally you shouldn't find that. It was a bit surprising when I found that even in a modern sample, because it meant that probably contamination was already happening in modern samples. When working with ancient samples, it was even more likely to happen, even though we take a lot of precautions, they're out there, the samples at first. Relying only on host DNA wasn't the only option. Also, here, they say that it's unlikely that the pig ate the human, but the opposite possibility is much more likely, humans eating pigs. If for pigs, it's not so much of an issue in the archaeological context, when archaeologists look at ancient feces, it has been much more often a problem to differentiate human from a dog poop. The problem is even doubled because in some civilization, it is known that people ate dogs. You would expect to find a mixture of human and dog DNA.

[9:34](https://youtu.be/gU4jx1pb8Tw&t=574)
The second step that we took to circumvent this issue is to look at the microbiome composition by using taxonomic profilers, such as Kraken2, then computing some sample pairwise distance metrics and doing some dimensionless reduction where it can get this nice plot. You can see your different samples in this dimensionally reduced space, in blue you have your dogs, in red and orange you have your humans, and in green you have some soil samples, you can see that they separate quite well. Based on this composition, and by comparing them to reference sample, here in this example, and say we have an imaginary unicorn gut microbiome profile, a dragon gut microbiome profile, you get this profile and you're asking, which profile does it look the most similar to? With some machine learning, a tool that's called SourcePredict, I realize I messed up the slide here. Okay, so sorry. Well, sorry, you're not going to see the slide I forgot, there... it's messed up. But with machine learning, you can identify the source of your sample by comparing them to reference samples. At the end... I messed it up so bad... At the end, you get a nice report that I can show you here. Okay. Sorry for that.

[11:32](https://youtu.be/gU4jx1pb8Tw&t=692)
Yes, so you get this interactive report, we have the summary table of the different findings. You have your microbiome embedding, so your samples in the dimensionally reduced space. Here are test samples, the sink samples are more or less within the human cluster. Then we get the summary plot where you have both the endogenous human versus dog DNA and the microbiome profile that are summarized in one single plot. We see that for two of the samples, it's quite clear that they are homo sapiens. For some of their sample either because there is a disagreement between the endogenous DNA competitive mapping and the microbiome profile or for some other reason, it is less clear. That's the summary report that you get at the end of the coproID pipeline, plus a lot of other things that I didn't mention that are specific to ancient DNA.

[12:45](https://youtu.be/gU4jx1pb8Tw&t=765)
The last release of nf-core/coproid was published in April 2020, so more than a year ago. More than a year in the nf-core history, short history, is quite a long time actually. That's why I put this picture here because if you look at the code of coproID, it looks like you're doing archeology of Nextflow code when you look at coproID. It was using nf-core tools version 1.8. In 1.8, there was no nf-core schema, and of course there was no nf-core DSL2 and no modules. When I look back at it, it's quite exciting to see that nf-core develops so rapidly, but it's also challenging to keep a pipeline up to date with the latest template and the latest evolution of nf-core, especially when you don't need to update it so much anymore because it's working, but you don't want... you don't need to add extra new functionalities. That's it from my end. The repository is nf-core/coproid. We published it in 2020, so the article is here, and there is a Slack channel, #coproid. If you have any questions, now is the time, and thank you very much again for inviting me to present coproID.

[14:33](https://youtu.be/gU4jx1pb8Tw&t=873)
(host) Thank you very much, Maxime. If you have any questions, you can unmute yourself, and so I've enabled that. You can unmute yourself and ask them directly. Alternatively, you can put them in the chat, and I can read them out. Let's wait a couple of seconds. Okay, so I don't see any questions pop up, and nobody has unmuted themselves. Thanks again, Maxime.

(speaker) You're welcome.

(host) We will be sharing the slides that Maxime presented today. They'll be uploaded to our website after being put up on Figshare. Now I'd like to announce that we have two more sessions lined up for you before the winter break. Please keep an eye out for announcements on the bytesize channel or on Twitter for future talks that will be coming up in January. See you next week.

</details>
