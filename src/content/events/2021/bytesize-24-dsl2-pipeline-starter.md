---
title: 'Bytesize 24: Where do I start writing my own DSL2 pipeline?!'
subtitle: Harshil Patel - Seqera Labs, Spain
type: talk
startDate: '2021-10-19'
startTime: '13:00+02:00'
endDate: '2021-10-19'
endTime: '13:30+02:00'
youtube_embed: https://youtu.be/Z_uPj7fAes8
locationURL:
  - https://youtu.be/Z_uPj7fAes8
  - https://www.bilibili.com/video/BV1nq4y197FR
  - https://doi.org/10.6084/m9.figshare.16836616.v1
---

# nf-core/bytesize

Join us for an episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 24: Where do you begin? DSL2 pipeline introduction

This week, Harshil Patel ([@drpatelh](http://github.com/drpatelh/)) will present an introduction to developing pipelines in Nextflow DSL2 using nf-core community standards.

Slides:

<div class="ratio ratio-16x9">
    <iframe src="https://widgets.figshare.com/articles/16836616/embed?show_title=1" width="568" height="351" allowfullscreen frameborder="0"></iframe>
</div>

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://youtu.be/Z_uPj7fAes8&t=1)
(host) [...] talk which is in preparation for the big DSL2 hackathon next week and today we have Harshil Patel from Seqera Labs presenting about writing your own nf-core DSL2 pipeline. During the talk if you have any questions for Harshil please put them in the chat and I'll read them out. As a reminder this video and all the previous bytesize talks are on the YouTube channel so you can consult them there later as well and all the links are on our homepage, so please take it away Harshil.

[0:26](https://youtu.be/Z_uPj7fAes8&t=26)
Hello everyone, good afternoon and thank you for watching and joining this talk. I will be giving the 24th nf-core bytesize talk. It's just amazing how we have got this far with this many talks now. I remember when we were still setting this stuff up, but hopefully they're useful to you guys. This particular one will be about trying to deal with the exploding head that is "I've got a DSL1 pipeline or I'm comfortable with DSL1, how do I now switch to writing a DSL2 Nextflow pipeline?".

[1:08](https://youtu.be/Z_uPj7fAes8&t=68)
Nextflow, as some of you know has this new syntax called DSL2 which is a more modular syntax. It allows you to reuse components and essentially make things a bit more flexible in terms of pipeline development as well. There's loads of resources so as I mentioned we have bytesize talks, that are revolved around DSL2 pipeline development. How we are tackling DSL2 on nf-core, but there are some generic concepts there, that may be useful to you as well. If you like, in the bytesize tests there are also some guidelines and other stuff that you may be able to take away to use for your own DSL2 pipeline development in non-nf-core pipelines too. Last week, Rike, Maxime and I recorded a couple of more up-to-date talks about the pipeline structure that we have for DSL2 in the template that we maintain on nf-core tools. I just briefly went through that again in preparation for the hackathon, so people get an idea as to what these files are doing and why they're there. Another talk which was typically a 45-minute talk. By my standards, that is where I typically get with these talks. That one was more about the process of picking a module and using the various tooling that we've created, to then contribute that back to nf-core modules. There are some generic concepts there that may be useful to you to take away.

[2:43](https://youtu.be/Z_uPj7fAes8&t=163)
The first thing I would do is attempt to watch some of these talks here, and also we've got loads of documentation and guidelines on the website, that'd be worth looking at. Figuring out how to tackle DSL2 in some way philosophically, would be in comparison to the way that you would attempt to figure out how these components are put together. The smallest unit of a DSL2 pipeline will be a module. So attempting to figure out how you need to mentally create and write these modules, to structure an entire pipeline together, is quite useful. We've got these guidelines and adding new modules and some of the concepts that we've used to attempt to standardize these these modules and the syntax available on the website for you to have a look at.

[3:33](https://youtu.be/Z_uPj7fAes8&t=213)
How you convert your pipeline will obviously depend on what status it is and where you're coming from, where you're starting from. For existing nf-core pipelines a lot of the groundwork has now been done, in terms of porting our old DSL1 template to a DSL2 template and switching things around, adding tests, all the boilerplates. This is one of the big advantages of working in a community like this, because it's not just me or another person, or someone else doing this, there's an entire community contributing to this. It just makes it easier to maintain and push these changes and across an entire sway the pipelines.

[4:15](https://youtu.be/Z_uPj7fAes8&t=213)
Those starting out with existing pipelines, you just would need to merge in the template sync that you would have got via a pull request. This is automated, it's sent out when whenever we release the nf-core tools package. You would have a pull request that is sitting there, waiting to be merged in. The first thing I would do with that, is just merge that in. You may have quite a few merge conflicts, but unfortunately, because this is such a big change, you will have to wrestle with those for now. Hopefully in the future they'll get smaller and smaller as things stabilize.

[4:45](https://youtu.be/Z_uPj7fAes8&t=285)
Get the tests working again. With nf-core pipelines we insist on having test data sets for continuous integration and for local testing. It just means that whenever you update the pipeline code you can test the pipeline to make sure it's working. After you've merged in the template sync, try and get the test working again, even if it's still in a DSL1 format. Just try and get the test working again and then decide how you want to tackle the implementation, which i'll go through briefly in some of the following slides. With new nf-core pipelines we have loads of guidelines and docs. There's very small differences in the way that the template has changed between DSL1 and DSL2, in order to manage the modules aspect. More change probably has been in the way that we've siloed away some of the boilerplate code into into lib directories, to make it easier to read the code that we've got there, and to update it. More guidelines and documentations on the website. I won't go through that in detail. Have a look.

[5:46](https://youtu.be/Z_uPj7fAes8&t=346)
But most importantly, if you want to contribute a pipeline to nf-core or you're thinking of doing it - whether you start off by with a new nf-core pipeline, or you start off thinking "I'm just going to write my own pipeline but maybe I might contribute to nf-core in the future" - please come and approach us first. Because we try not to have redundancy in pipelines and it's always a community decision as to what gets in. It'd be great if you can approach this before you lay down in your code. Things can become a bit awkward, when you've written a pipeline and it's fully functional and it's in all its awesomeness, but it may not necessarily fit with what we have or what we require. In that case we don't want to disappoint you or make things difficult for everyone.

[6:34](https://youtu.be/Z_uPj7fAes8&t=394)
If you've got a non-nf-core pipeline, so for example you don't want to contribute to nf-core, you can still use the tooling that we've created. It's completely up to you and it's flexible, how you adopt the standards and the template and the different files and stuff that we've got in the template. This is simply done by using the `nf-core create` command in the nf-core tools package. It's one command, you get a bunch of boilerplate stuff, that you don't have to do yourself. If you're just looking at writing simple Nextflow pipelines this may be overkill. But if you're seriously thinking about writing your own pipeline, even to use as a reference to see how the community itself is adopting best practices, how they run github actions, use continuous integration, configuration, linting, all sorts of other stuff as well. This is very useful to have a look. Like I mentioned, I gave a talk about that last week, so you can see what it looks like and what the files in that repository are doing.

[7:28](https://youtu.be/Z_uPj7fAes8&t=448)
It also means that you can sync in the template. Whenever we do a release of nf-core tools all nf-core pipelines automatically get this sync PR. Anything we've updated in the pipeline template then gets pushed automatically to these pipelines via this sync PR. If you have created a pipeline template, or be it, whether you're not going to contribute to nf-core, then you can use the sync functionality to update your template too with that. It just it just allows you to keep up to date with the best practices and other boilerplate and bugs fixes and stuff, that the community is implementing.

[8:07](https://youtu.be/Z_uPj7fAes8&t=487)
It also means the pipeline can be contributed later. Like I said, approach us first if you are seriously thinking about it. When you use `nf-core create`, it does a few things, especially with git. There's a bit of magic there, that allows you to then contribute that pipeline to nf-core later on down the line, if you so wish to so. That's another advantage. We also have loads of other nf-core tools commands. Check it out. I won't go through them now but there's various tools for linting and other stuff that'll be useful for maintaining and developing the pipelines.

[8:44](https://youtu.be/Z_uPj7fAes8&t=524)
The first call I would recommend is to look at nf-core modules. That's our repository for wrapper scripts, essentially, or DSL2 modules. It's been developing immensely well. We've got six to seven contributors. We've almost got to 300 modules now, which after the hackathon, I imagine, will completely surpass that. It's just a repository for standardized module wrapper scripts for individual tools like fastqc, or trim galore, that you can just pull and use directly in your pipelines. You don't need to go through the effort of writing these modules. It saves a lot of work and this fits in with the ethos of Nextflow DSL2 as well. It's constantly evolving. I won't say that it's completely stable, because it's not. I would say however that we're constantly making it better and trying to shift towards using as Nextflow-esque language and approaches as possible. So check it out!

[9:47](https://youtu.be/Z_uPj7fAes8&t=587)
To add to that, we've got loads of tools that we've added in nf-core tools, specifically to deal with modules. I've listed them here. You can list modules, install them, update, and all sorts of other functionalities. Some of this I refer to in the other pre-hackathon talk about contributing to nf-core modules, that I gave last week. The link was in the first slide.

[10:10](https://youtu.be/Z_uPj7fAes8&t=610)
We plan to have subworkflows in the future. For those of you that don't know what a subworkflow is, it is essentially a chain of modules. A module is a unit of DSL2, let's say where you've got FastQC that runs on a single sample and performs a particular task. However you can chain these together so you can run FastQC and adapter trimming after that as a subworkflow. Then you get a larger chain of modules, that you can then just plug into a pipeline, without having to individually chain them together. This is the true power of DSL2.

[10:48](https://youtu.be/Z_uPj7fAes8&t=648)
How we make subworkflow shareable and reusable across pipelines is going to be a real test, because we'll have to figure out a few other things. We plan to tackle some of this at the hackathon next week and maybe - as you can see we've got an nf-core modules command - we'll probably have an nf-core subworkflows command, that will install all of the module dependencies, as well as the subworkflow, wherever it needs to be installed. All you really have to do is include it in your pipeline.

[11:16](https://youtu.be/Z_uPj7fAes8&t=676)
Getting started, I would probably start by looking at existing pipelines that have done this. I'm pointing to nf-core pipelines here, because it's what we know, it's what we've done. There's a full list. If you click on that link there, there's a DSL2 tab on the pipeline health page on the website, that allows you to look at other example pipelines, if they're more applicable to you. I think most importantly is setting up a nice test data set. We try and tackle that right at the beginning. It's always good to test your pipeline right from the offset. It also means that other people can collaborate on the pipeline with you and you can identify bugs and issues and pull requests or locally, that you can fix, whilst developing the pipeline together.

[12:05](https://youtu.be/Z_uPj7fAes8&t=725)
I would say it's incredibly vital to have a nice minimal test data set that you can use. And also this becomes important when your people just want to test the pipeline on their own infrastructures, for example. This minimal test data set is independent of the samples they're using and so you know the test data sets should be working and it allows you to rule out other issues with infrastructure and such, when using Nextflow.

[12:30](https://youtu.be/Z_uPj7fAes8&t=750)
Compiling list of modules. This is quite an obvious point, but you need to know what modules you have. We've got loads on nf-core modules, like I said we've got already got about 300 of them. A lot of this work has probably already been done for you. That's not to say, we wouldn't like your contributions there too, because then it just means it's done for someone else as well. Hopefully at some point we'll get to a point where... a majority at the moment is quite genomics focused, but hopefully we'll be getting other modules in there, from other life science areas as well.

[12:30](https://youtu.be/Z_uPj7fAes8&t=750)
Find and recycle sub-workflows. We're still working on this, or adding this to nf-core modules, or having a separate repository, maybe nf-core sub-workflows for this. These are at the moment within pipelines like rnaseq. The sub-workflows folder. You can have a look in and see if there's anything you'd like to reuse from there or from anywhere else. At the moment it's a manual process but hopefully we will automate this in the future. How you do this and collaborate on these modules depends entirely on how you want to develop the pipeline. You could create a list of modules as separate issues and then work your way through those or you could create a project board, like Sarek has done, which I will show you in the next slide, or maybe the side after. You can collaborate and tick these off the list, eventually, whilst you're developing the pipeline.

[13:54](https://youtu.be/Z_uPj7fAes8&t=834)
In terms of implementation, we've built a lot of these tools. Nf-core `modules create` is an is an example of this. It just takes a vanilla module template with loads of to-do statements and other things in, that are really useful for newbies and beginners. Just as a reminder, to make sure that you filled in the correct bits in the file. It has a load of to-do statements within this particular template. When you run `nf-core modules create` it just replaces the name of the module that you'd like to create within this template. Then you have to go about then replacing the bits you want in order to finesse and add your module or create a module. Some of this stuff I went through in that "contributing to nf-core modules" talk, so please do have a look at that and you'll get an idea as to how that can that can work. You can do that for both local modules, the ones you don't want to contribute to, or ones that you do.

[14:52](https://youtu.be/Z_uPj7fAes8&t=892)
Reuse biocontainers. The biocontainers are essentially bioconda packages built within both Singularity and Docker containers. It's an awesome resource, that we've been using almost exclusively for all of our modules. It just means that you can get a Docker container and a singularity container for free. We don't have to maintain anything. If someone adds a new bioconda package, we get that as a container for free. Reusing this is nice because then it just gives you this option and of not having to host and maintain this yourselves. Passing sample information around is also quite important. You need to figure out the flow of your pipeline. Typically what you would do is have you'd have different values in a channel for different sample attributes, but this gets a bit complicated when you want to generalize a module. The best way to do that is to put all of this sample information into, what we have called, a meta map. Then you can have as many - it's like a python dictionary - you can have as many attributes within there and pass that through a pipeline. That also means you can reuse existing modules and nf-core modules and so on, and still have access to that meta within your pipeline context. How you do that use is also something you need to think about.

[16:06](https://youtu.be/Z_uPj7fAes8&t=966)
Try and stick to a single syntax convention. We obviously have our own, but whatever you do, whether you want to develop your own, just stick to a single syntax convention. Because it just makes things consistently easier to maintain and to update over time. When you change that syntax, or that convention, you can reuse what we've done. A lot of people are. The caveat there is that it's constantly evolving, it's something you have to keep on top of. I don't think that's a bad thing, personally, because everything changes, everything evolves. Please write your modules in a way that they can be reused! That's the true power of all of this. That just means that someone,whether you contribute them to nf-core modules - whether you write them and keep them locally within your pipeline - it means other people can just pull that module straight away and reuse it without having to do much.

[17:00](https://youtu.be/Z_uPj7fAes8&t=1020)
There are various different approaches in terms of how you tackle the implementation. I particularly prefer the bottom-up approach, where you have your main script. It's just completely DSL1, you essentially just comment out the whole thing and start adding one-by-one each of these modules into the pipeline. It just allows you to test every step of the way. Also there's other things around the way - that you pass channels and manipulate channels to these modules, that it allows you to do quite interactively, whilst you're developing the pipeline. This is what what I prefer doing. There are other approaches. For example, there's a couple of links to issues there, where we've been creating a list of modules, that you can see on the right here and that we have worked our way through. This can be individual issues or you can create a project board which is what that nf-core/sarek link will take you to.

[17:54](https://youtu.be/Z_uPj7fAes8&t=1074)
You can do a top-down approach, where you write your modules first and then stitch the pipeline together. There may be caveats in the way that you do that. It's not impossible. I think Praveen and Maxime have done that with nf-core/rnavar, where they wrote the modules and then stitch the pipelines together. There are a couple of caveats in terms of the way that you may need to update modules, whilst you're then developing the pipeline, because you've already written them and you may need to change them, to fit into the pipeline and stuff. It's still a valid and plausible approach.

[18:28](https://youtu.be/Z_uPj7fAes8&t=1106)
We will be changing the syntax that we're using for DSL2 very soon, hopefully. We're moving to a more Nextflow native syntax. I've provided a brief description of this in the "contributors to modules" talk at that particular time, if you want to skip through to it, so I won't go through this in any detail. The information is there and hopefully it will make the adoption and usage of these modules even more widely accessible, because we'd be using a native Nextflow syntax. That just removes things like the functions file and other things, that have been a bit of an issue in terms of customization. Watch this space. Things are evolving. Everything will be updated and hopefully you'll be able to keep on top. This is something that we'll probably be discussing and trying to iron out at the hackathon.

[19:23](https://youtu.be/Z_uPj7fAes8&t=1162)
If you need to get in touch: Slack. There's the #modules channel on there, we have a DSL2 pipelines channel somewhere as well, I can't remember what it's called. James will probably tell you when you finish. We've got we've got another channel, for those that are... DSL2 conversion channel or something like that, for people that are interested in knowing more about the conversion process. Github, Twitter, Youtube. Reach out however. Join the community, join the slack workspace! There's a lot of information to be gained there and it's ridiculously easy to join. Even if you're just lurking in the background you'll be surprised how you can just absorb the knowledge. Thank you to everyone in both communities, nf-core, Nextflow, bioconda, biocontainers, my new worklord Seqera labs. An awesome bunch of people and that I will be seeing in Majorca tomorrow. We have a hackathon next on the 27th - 29th, which has partly been mentioned a few times now. If you haven't signed up, I think the sign-up is still open. I look forward to seeing you there. Thank you.

[20:31](https://youtu.be/Z_uPj7fAes8&t=1231)
(host) Thank you very much Harshil. The channel you were talking to talking about is slack channel is #dsl2-transition.

(speaker) There you go. That's it.

(host) Are there any questions? You need to post them on Zoom or in the Slack chat. Anything? Well, I guess that's it for today. Of course we have the hackathon next week. That will be in Gathertown, so it'd be a very nice interactive environment. I'm sure you can come by and ask Harshil questions then. But otherwise we do have a bytesize talk next week from Daniel Straub talking about nf-core/ampliseq. This will be normal time: one o'clock CET on tuesdays, and then like Harshil said, we have the hackathon Wednesday to Friday next week.

</details>
