---
title: 'Bytesize: resources to learn Nextflow'
subtitle: Sateesh Peri - Bioinformatics Pipeline Engineer, Bluestar Genomics
type: talk
startDate: '2022-05-31'
startTime: '13:00+02:00'
endDate: '2022-05-31'
endTime: '13:30+02:00'
youtube_embed: https://www.youtube.com/watch?v=oO7rAp-QZOk
locationURL:
  - https://www.youtube.com/watch?v=oO7rAp-QZOk
  - https://doi.org/10.6084/m9.figshare.19960847.v1
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: resources to learn Nextflow

This week, Sateesh Peri ([@sateeshperi](https://github.com/sateeshperi)) will talk about recources that are available for newcomers to Nextflow - but even more experienced Nextflow users might learn something!

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=oO7rAp-QZOkt=1)
Hello, everyone. My name is Franziska Bonath. I'm today's host of this talk, and with me is Sateesh Peri. And he is going to give a talk about resources to learn Nextflow. In the chat, he put a link for live slides. And from now, I give it over to Sateesh.

[0:26](https://www.youtube.com/watch?v=oO7rAp-QZOkt=26)
Thank you, Franziska. As mentioned, we have pasted a link to the live presentation that you can follow along. I'm going to be sharing my screen as well. So today here, I'm going to give you a walkthrough on all the resources that are currently available to get a good handle on the Nextflow workflow management system.

[1:01](https://www.youtube.com/watch?v=oO7rAp-QZOkt=61)
Lets dive straight into it. If you're here, some of the goals that you should be having are to learn a simple syntax for writing pipelines that enable you to reuse existing scripts and tools from past prototyping, and also to develop self-contained pipelines, manage versions, and to be able to reproduce any form of configuration on demand. One of the most useful resources was this blog post earlier in the year, that was very useful in bringing together all the useful links that are available. But since then, there have been a lot of developments, DSL2 is now default. We were wondering if we could give a more curated list of resources that are available to learn Nextflow through self learning. Before you dive into Nextflow and nf-core, some of the prerequisites that you need to master are scripting languages, bash programming, some experience with containers and software dependencies, and also quite some experience with version control systems. Not only for version controlling but also collaborating with open source communities such as nf-core.

[2:34](https://www.youtube.com/watch?v=oO7rAp-QZOkt=154)
Here are some links. There are enough links online to help you get up to speed. But here are some links that I have followed. Pretty much this presentation is going to be about how I have learned Nextflow over the past year. I'll be sharing links that I have found useful throughout my journey. These are some of the links here to get you started on all the prerequisites that you might require as you are about to get started learning Nextflow. Once you have those down, well, you're here because you want to learn Nextflow. Nextflow enables scalable and reproducible scientific workflows using software containers. It allows the adaptation of pipelines written in most common scripting languages, and its a fluent domain specific language, currently in its version two of evolution, simplifies implementation and deployment of complex parallel and reactive workflows on cloud and hybrid environments.

[3:45](https://www.youtube.com/watch?v=oO7rAp-QZOkt=225)
The core features of Nextflow that should interest bioinformations is it enables workflow portability and reproducibility, simplify parallelization and large scale deployment, easily integrates existing tools, systems and industry standards. These are just some of the features. I have here linked two video links that give an overview of nf-core and Nextflow and also some of the details of why nf-core or Nextflow stands out in comparison to other workflow management systems. I would recommend if you're a complete beginner to go over these videos as well. But once you're through, we are talking about you might already have an introduction to Nextflow. But I want to like completely go through the nuts and bolts of to learn the whole syntax.

[4:53](https://www.youtube.com/watch?v=oO7rAp-QZOkt=)
We are recommending thrse courses for you here. Because now DSL2 syntax is default, we would like to transition most of the courses to DSL2 syntax as well. In that aspect, we are recommending three courses here. One is the Nextflow tutorial variant calling edition that was developed by myself and colleagues at CDC. This particular workshop content was developed so that it explains the concepts through a variant calling workflow example. I'll briefly walk you through the contents on how it has been divided for this particular workshop that we have developed. It starts off with an introduction to Nextflow, and then we dive right into introduction to nf-core with the main aspect being that we wanted to tell participants not to reinvent the wheel in writing pipelines. And an introduction to nf-core is also a great way of introducing how to run Nextflow pipelines straight out of the box, and also to show how they can be done on local & on-prem HPC clusters. We introduce to these sessions here, and then we dive into the details or start with the syntax of Nextflow, then get into the details of channels, processes, workflows, and operators. And finally, we give a challenge of converting a bash script of a variant calling workflow into Nextflow processes. And finally, we show the workflow. Then content will progress into modularizing all of the processes that are involved in the variant calling workflow and towards the end is where we have introduced a chapter on assembling the variant calling workflow now using nf-core modules. Previously we show how you can develop your own local modules, but we also show how the same process can be done using nf-core modules. It comes as a full circle towards the end. This session, we have designed for nine hours and it has worked out pretty well in training.

[7:43](https://www.youtube.com/watch?v=oO7rAp-QZOkt=463)
The next course that we would recommend is the Software Carpentry RNA-Seq workshop content that is currently still being developed. It's in the pre-alpha stages, but it's a great resource as it has more exercises and also it deals with an RNA-Seq example versus the variant calling example that we have previously dealt with. Further we have the official training material from Seqera labs. This is much more organized and directly from Seqera labs, but they're still in transition to move from the DSL1 examples, as employed in the current version, to DSL2. There'll be updates soon to look forward here.

[8:34](https://www.youtube.com/watch?v=oO7rAp-QZOkt=514)
For all these courses, we are trying to offer students to use Gitpod environments. Gitpod is an open-source developer platform that quickly spins up a virtual environment from a Git repository. And you can have it pre-installed with software such as Nextflow, Conda, and Docker. All the courses that I have just shown you, at least the variant calling workflow and the Seqera labs training material, can be opened in Gitpod. The links are accessible in the setup pages and the workshop can be followed in the Gitpod environment. Spinning up a Gitpod environment and walking through the tutorials at your own pace should give you a head start.

[9:37](https://www.youtube.com/watch?v=oO7rAp-QZOkt=577)
The next things that you should be bookmarking for handy access are the documentation links. These are the pages that are definitive in terms of the changes or the processes associated with Nextflow or nf-core. These should be your go-to docs at any time for anything definitive about Nextflow and nf-core. Bookmarking these is a thing. At this point, you might be at the stage where you say: I know Nextflow and I'm ready to take on the fight. This is where your next point is to introduce yourself to the community, the nf-core community, where you'll find a lot more talented and expert people. Nf-core is this community that's using the Nextflow workflow language to build a set of curated, peer-reviewed best practice pipelines. Nf-core pipelines have strict guidelines. If one of them works, all of them will. And all these features make nf-core pipelines, even if you're not a developer, if you just want to use the nf-core pipeline straight out of the box, these features absolutely make it possible.

[11:12](https://www.youtube.com/watch?v=oO7rAp-QZOkt=672)
In your next step of journey, once you have finished the courses, a lot of time we would recommend spending and looking over at nf-core available pipelines, the nf-core modules that have recently been integrated into the web page itself. And also to look at the nf-core tools that are available if you're interested to develop your own pipeline in the future. Throughout, as you're going, also a great resource just like this bytesize talk, are our other bytesize talks pertaining to any specific topic. There is a YouTube list. As you go through different course materials, you might want to check out those specific videos associated. You want to make sure that you are looking at videos with the transition to DSL2 in mind. Some of them do still have DSL1 syntax. Just so that you know.

[12:24](https://www.youtube.com/watch?v=oO7rAp-QZOkt=744)
The next great resource that you should be spending a lot of time, especially getting in touch with the community is through Nextflow and nf-core Slack channels. Now, these are an absolutely great way of connecting with people around the world and most importantly, learning from them. You will see in the Slack channels, we have different channels pertaining to different topics. One pertaining to modules and request review, just as an example here, where you can actually ask questions. And more importantly, for me, Slack has become a way of troubleshooting/debugging. Any eroor messages that I have, the first thing is I'm going and checking in the Slack workspace by searching for it. And most of the time somebody has already asked for it and there's solutions already recommended from the community. Slack is a great way of learning a lot of things and also looking for help in troubleshooting. The Nextflow help channel is a great place for asking help from everyone in the community. It's great that there are people from all over the world, all over the time zones that can help you at all times. Definitely make use of this resource. And if at all you are in a capacity by the end where you have mastered the skills yourself, please give back to the community. There are always newcomers coming in and it's always great to have their questions answered by somebody else in the community.

[14:18](https://www.youtube.com/watch?v=oO7rAp-QZOkt=858)
The next learning resource is GitHub. The GitHub review process, it not only makes sure that the community coding standards are adhered to, but also it's a great way of, again, learning the whole process. As you have experts from the nf-core core team guiding you on making sure all the standards are met too and also giving you unique solutions. Interaction through GitHub is an important part. As you are interacting with the nf-core community, just pouring over the existing GitHub repos of nf-core pipelines and looking at the processes are great places to learn as you start building your own pipelines. There's more to collaboration than you think. One of the resources that I would recommend is the Turing Way guide for collaboration. Especially if you are new to GitHub and would like to know more about how the review process works and all that you can find more information through this guide. At this point, you should be quite capable of doing any of these things especially.

[15:49](https://www.youtube.com/watch?v=oO7rAp-QZOkt=949)
One is to run the nf-core pipelines straight out of the box. Especially if you're a facility and there's an nf-core pipeline that fits your needs. One is to use the nf-core pipeline straight out of the box. Second is to run the nf-core pipelines with some modifications. Either you turn off some of the features or you add some features or modules on your own with some modifications, or create an entirely new pipeline from scratch using the nf-core template. Now using the nf-core template gives the advantage of making sure that they all have the same structure as the other nf-core pipelines. You can then adapt a lot of pieces, especially modules and subworkflows from other nf-core pipelines into your templated pipeline that you're creating. In that case, all you have to do will be to basically shop for modules in nf-core if you have all the modules available there or consider making one if it isn't there, so that you are giving back to the community as well. And in some cases, you'll have to keep some local modules.

[17:10](https://www.youtube.com/watch?v=oO7rAp-QZOkt=1030)
Once you plan all this, you can just build your pipeline from the ground up.

"Through open discussion and collaboration among the community, it's possible to leverage the knowledge of experts across the world for the development of domain-specific pipelines and implementation of current best practice analysis methods."

This is from the nf-core Nature paper, and I truly think these words hold value and I have seen it in practice. I have learned from this community a lot, so it is possible and we should keep doing this. And just so that you know, there's a new mentorships program being launched jointly by the Nextflow and nf-core community. They're hoping to organize mentor and mentee pairs, and especially reach underrepresented groups and areas to improve outreach. You'll find more details about this on the nf-core website. I think the current round is already done but they have a couple more rounds coming as well. With that, I just want to say thank you and welcome to the community. This is pretty much my setup that you can see here but I'm sure many of you will resonate with that as well. This is currently my short guide on resources to learn Nextflow. We'll soon have more improvements, but we'll be back with those again, thank you so much everyone.

[18:56](https://www.youtube.com/watch?v=oO7rAp-QZOkt=1136)
(host) Thank you very much. I now open the floor for any questions. You can either give a hand sign or just come in. There doesn't seem to be any questions, then I would like to thank you again and for anyone who might have a question, you can always go to the Slack channel as we just discussed, and ask your question there. And also I would like to thank the Chan Zuckerberg Initiative for funding of this, and I hope to see you all soon in Slack and GitHub and everywhere else. Thank you very much.

(speaker) Thank you all.

</details>
