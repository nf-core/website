---
title: 'Bytesize 5: DSL2 module development'
subtitle: Harshil Patel - Francis Crick Institiute, UK
type: talk
startDate: '2021-03-02'
startTime: '13:00+01:00'
endDate: '2021-03-02'
endTime: '13:30+01:00'
youtubeEmbed: https://youtu.be/ggGGhTMgyHI
locations:
  - name: Online
    links:
      - https://doi.org/10.6084/m9.figshare.14160116.v1
      - https://youtu.be/ggGGhTMgyHI
      - https://www.bilibili.com/video/BV1aK4y1D7z3
---

This week, Harshil Patel ([@drpatelh](http://github.com/drpatelh/)) will present: _**DSL2 module development.**_ This will cover:

- Module file structure
- Writing new modules
- Automated testing

The talk will be live-streamed on YouTube:

- YouTube: <https://youtu.be/ggGGhTMgyHI>

<details markdown="1"><summary>Video transcription</summary>

:::note
This text has been edited to make it more suitable for reading.
:::

[00:00](https://www.youtube.com/watch?v=ggGGhTMgyHI)
Hi everyone, so this is our fifth byte-sized talk which I think have been really really useful of late, especially sort of summarizing what we're doing on nf core in small bite-sized chunks so people can get more familiar with what we're doing. It also acts as sort of a persistent archive of how to do things on nf-core and and hopefully also with Nextflow and git like Alex's talk [https://www.youtube.com/watch?v=gTEXDXWf4hE] last week and other things so thank you all for joining.

[00:33](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=33)
Today i'm going to be talking to you about nf-core/modules which is our Nextflow DSL2 effort so, as some of you may know, Nextflow now (actually last year July I believe it was) released a new modular syntax called DSL2 and that allowed a lot of flexibility in terms of pipeline development and also it got us very interested here on nf-core because we have 40-50 pipelines that share functionality and do similar things and we try and standardize this as much as possible to help with the development of these pipelines and also other users to for them to be able to understand what they're doing and things like parameter names and configs and all of that sort of stuff. So, this was actually really a really big thing for us.

[01:27](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=87)
To take you straight to some sort of terminology – What is a module? So, a module in our definition is something that is as atomic as possible it can't be broken down into anything smaller. So, you would imagine fastqc being an example of that as it’s a single tool to perform a particular task and that would then be termed what we would call a module. Similarly, you might have BWA-MEM or BWA Index; these are all single tools that perform a particular task and in our definition that is what we would call a module.

[02:01](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=121)
You can also have sub workflows which are chains of modules that you can use to offer some sort of functionality within the bigger term of the word which is the workflow and a sub workflow would be something like sorting a bam file, indexing it, then running some stats on it and having all of that packaged up into one sub workflow as a chain of tasks and you can imagine, in genomics especially, you would use that sort of sub workflow quite often because you would create multiple iterations of a BAM file doing filtering and marking duplicates and doing other stuff and at each of those points you may want to sort index and run some stats on it. So, actually I think the most powerful aspect of DSL2 will be sub workflows, well written sub workflows.

[03:19](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=199)
A workflow is an end-to-end pipeline, so in DSL1 language it would be a pipeline that runs from end to end and back with DSL1 there was a lot less flexibility in terms of what you could do with the workflow and how you can include things, how you can overwrite parameters and so on whereas now with DSL2 this has become a lot more flexible.

[03:50](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=230)
We can share the DSL2 modules between pipelines if they're well written well enough and we can even share sub workflows between pipelines. As I mentioned, this is sort of done manually at the moment but, hopefully, in the future we will have some functionality to deal with that.

So, when we went about trying to figure out how we would deal with nf-core/modules and it took us a while of procrastinating because it's not a trivial task making wrappers like this as standardized as possible for an entire community of people to use and also for them to be flexible enough so you know you're not imposing certain settings and options onto developers themselves (You may want to use fastqc with different options or want to publish it in a different directory and so on).

[04:37](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=277)
So, a number of key things that we don't compromise on nf-core, if you know nf-core - Reproducibility was one of the earliest things we discussed, how would we make these things reproducible. At the moment, the way it stands now is that the module itself is installed physically within the pipeline repository so when you do a git release or github release, that module will be shipped with the pipeline. So, in that aspect, you can always ensure that the module is reproducible. You want moving that module anyway, you won't do anything, you've got a static representation of that module within the release code and so that's where the reproducibility comes from.

[05:16](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=316)
We're also more recently decided that we would use Biocontainers for all of our software packaging. We initially started off by using or building docker containers for using environment YAMLs and so on, like similar to what we're doing with DSL1 and the nf-core pipelines, but in the end, we decided that reusing Biocontainers is much more advantageous. We don't have to have an infrastructure to deal with that, we don't have to build Docker containers and the great thing now also is that with recent updates in one of the Nextflow edge releases you can also directly download Singularity containers. So, you don't need to convert the Docker to Singularity containers, which again, is another thing that we've traditionally been doing with DSL1 pipelines.

We have one Docker container that gets downloaded, converted to Singularity and that is what is used then by the pipeline but now, Biocontainers are also hosting singularity images directly so we don't have to convert anything. We're directly downloading them over HTTPS and using those.

[[06:20](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=380)]
And that's, that's amazing, because we've had a number of issues with users running out of whole space in home directories and so on and this sort of bypasses all of those issues. And obviously, supporting conda as well, which is where I guess you would imagine fundamentally these Biocontainers are built from. Biocontainers are essentially conda packages built in containers, so either Docker or Singularity.

[06:42](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=402)
You can also upgrade and downgrade these modules, if required, within the pipeline repository so this sort of imposes the restriction that they have to be relatively self-contained. When you install a module, which means that you can have different versions, um sorry not different versions, you either have a version of samtools using version 1.1.1.10 or you can have one using 1.11. It's completely up to you how you manage that and, in order to fulfill that sortof criteria, we need to have them be as flexible as possible.

[07:16](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=436)
Documentation was number one on the issue list for nf-core modules when I
initially created this back in, sort of, July 2019 and, again, documentation is quite key to all nf-core pipelines and it took a bit of thought but we've now decided on having, sort of, a yaml file that gives a brief description with tool input-outputs and and the authors that have contributed to it.

[07:40](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=460)
For simplicity and also for learning curve we wanted to stick to using the Nextflow coding style or the coding pattern that is familiar to most people writing Nextflow workflows in order to make it easier for them to contribute to not only nf-core modules but also to install the modules themselves and to figure out what's going on.

I think that's quite important so the simplicity there is incredibly important for the learning curve when it comes to figuring out what these modules are doing. I mean, me personally, I find DSL1 modules are really great, DSL1 pipelines are really great and that you have everything in one workflow because it's more findable that way. With DSL2, you can package things up and put them in various different places and it's not always that trivial to find them, so the way that we standardize the structure and and the way that where we're writing these Nextflow imports and so on is actually quite important.

[08:35](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=515)
Also, we're using certain parameters and other things and just generally standardizing how we're doing things across nf-core pipelines and, hopefully, this can also be reused by the Nextflow community. And the great thing also is, if you update it on nf-core/modules, where if you update the version of samtools on nf-core/modules because new releases come out, then everyone benefits from that so it becomes a bit like the way the Conda operates in updating builds of their packages and so on.

[09:04](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=544)
So, I've mentioned the standalone and why it's so important to have this in Nextflow workflows because you can just install any given module you want and it works. I mean, great thing about Nextflow also is the fact that it's built on top of Groovy which is it's own programming language and so we can exploit and use Groovy syntax, Groovy code and I've actually been getting into it a lot quite recently and some of you will notice in the release after next, where we, hopefully, will have a DSL2 pipeline template release where we've siloed a lot of the boilerplate code away into Groovy lib functions and so on and it just tidies the code up so much more.

[09:46](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=586)
But the advantage here is also that because Nextflow is built on top of Groovy, we can use Groovy to write functions and other things that we can use to manipulate or change things that are not as possibly as trivial with Nextflow.

[10:03](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=603)
So, automation, again, quite important, you know, we want a situation where pull requests submitted to the nf-core/modules repository and everyone that's reviewing that is also happy with the fact that whenever that pull request is created, we're running the right tests, linting tests and also, now we are running tests because these wrappers are self-contained, these modules are self-contained.

[10:29](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=629)
We have Conda, Singularity and Docker definitions physically within the main script of of the module and so we can now test using CI whether that module works with Conda, whether it works with Docker and whether it works with Singularity and as a layer on top of that, we've added this ability to use pytest workflow (which mainly was done by Edmund Miller) where you can actually generate the outputs from the module, create md5 sums and now, you're not only testing that the module works, you're testing that the module is producing the same md5 sums and this is quite important because it would be quite easily overlooked whether something is being changed through releases of a module and so on. So this becomes a really nice way of unit testing these modules and it's working quite well now and we're still sort of early days in it, but it's working really well.

[11:26](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=686)
The portability is mostly taken care of by Nextflow (Nextfow is amazing as you know), it works on virtually any platform. But there is, I think, a caveat here in that these modules have to be written as simplistic as possible to adhere to Nextflow guidelines on running on multiple platforms. We can't add customization and that would violate that essentially and so we've taken a lot of care in making these modules as simple as possible so in fact that they are portable on these different platforms.

[11:58](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=718)
And, of course, I made up a buzzword there because I couldn't think of one, um, but the last but not least, people have to be able to install these modules and use them themselves so they have to understand what they're doing and so again this comes back to the simplicity concept; as well as newcomers that may not have that much knowledge with Nextflow or nf-core pipelines for that matter and want to contribute to nf-core/modules. We've attempted to try and make that as simple as possible.

[12:25](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=745)
We had a number of ideas as to how we would deal with these modules, how would we physically use these modules within pipelines. A number of ideas were put around
back a year or two ago, where the folks from Bioconda got involved, a bunch of us got involved from here on an issue on nf-core/modules and we discussed the idea of using conda to manage all of these modules, get sub module, npm; but in the end we decided to go for something a bit more simplistic, which was using our existing nf-core/tools package that we already have to maintain, um you know, things like ’create pipelines’, to lint pipelines and so on.

[13:07](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=787)
We thought we'd add an extension to that which would be the ’nf-core modules’ command and that would do a lot of this stuff and in an overview, in a nutshell, it's actually a very simplistic approach, so, you know when you're installing a module all you're doing really is querying the github api and installing the module. We can now also add on other tools to allow us to lint to make sure that the standards that we've set for those modules in terms of syntax and other things that you normally might miss on a pull request even though, you know, you may be reviewing just a few files it's quite easy to miss that documentation has not been added in the right places and so on. So, we're in the process of extending that.

[13:45](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=825)
Some of this was already available, so the install and list commands are already available in in 1.10.2 and actually that's all I needed for the latest RNAseq DSL2 implementation but now we're adding a bunch of other stuff on top which is cool and a lot of this has actually been done by Kevin Menden.

[14:04](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=844)
Just to give you again, going a bit backwards, um, so the repos created 2019 just literally whilst I was at ISMB in Switzerland and, um, it was a moment of inspiration. And it kind of sat there for a while we didn't do much about it. I had a brief skeleton that I added. And then we’d done some at the Crick hackathon which I organized in March 2020 where Phil, Felix, Anna, a bunch of others started adding a few modules. Phil sorted out some Docker pushes for these modules and we made a bit of headway there. But I think the real dent we made was in the July hackathon, which was organized by Gisela, Enrique and that was our first remote only hackathon.

[14:49](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=889)
We really sort of sat down and stripped this apart and I think that's exactly what we needed, we stripped the modules into different components, had a lot of discussion as to how we would organize things how we would pass options around and, you know, reproducibility, which containers we would use. We made quite a big dent so after that I was relatively happy with the progress we'd made and we had a plan.

And, of course, a week after that Paolo released version 20.0.7.1, which was the first Nextflow release where you're enabling DSL2 and not previewing it. So, that meant we really had to do something about the modules and things that sort of made sense to me in terms of the way we attack this is that we have a proof of principle implementation in a real pipeline to see how it would work.

[15:46](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=946)
And, so, I went about step-by-step updating nf-core/modules and also completely rewriting the nf-core/rnaseq pipeline from scratch with DSL2 and so we released v2 of that pipeline and then there were some other issues with the alignments we were using and then the methods we were using, the pipeline even. And so, we done another release, v3 quite soon after that. But what that allowed us to do is actually gauge how flexible it would be for these modules to be used in a real pipeline, in a real end-to-end pipeline and that really sort of triggered a lot of this stuff and then after the rnaseq release, I then went and updated the nf-core/tools pipeline template with that, so eventually that will be released as I mentioned the week after, not the week after, the release after next, in the nf-core/tools as a template. Hopefully, all nf-core pipelines switch to that in the future.

[16:45](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1005)
Over the past month or two it's just been crazy, I've been working on re-releasing and converting the viralrecon pipeline to DSL2 and a number of people have stepped in and helped out with that, so thank you Jose, Kevin, um a bunch of others, Michael, Anders. And so, they've added modules so the repository sort of become bloated. A lot of these modules have been added in the last month or so I would say and in that process, we've also been refining the CI tests.

[17:14](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1034)
Kevin has also currently got a PR coming for using a very nice standardized minimal test data set. So, it just means that we can reuse the test data as much as possible across different modules without having to add the. And the more standardized we can make this the better it's going to be without having a thousand randomly named files in a repository that we're using for test data.

[17:37](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1057)
And also, as I mentioned earlier then, Kevin's been adding linting functionality, md5sum checks and there's a bunch of others that we're planning to add some point soon.

So this is what a typical module will look like on nf-core/modules. You have (get a bit fancy) the nf-core/modules top level directory, then you have a software directory, a module which in this case is fastqc, this functions.nf file which we're using to bring in some custom Groovy functions to deal with a few things in the main script and these are shipped with each module, so you have one per module.

You have the main script which is doing the crux of the work. It's just an excellent process, a single Nextflow process and you have a meta YAML as I mentioned earlier that documents them.

[18:29](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1109)
A brief description of what fastqc is, what the inputs are, the formats, the file extensions and also the author list. You also have a tests directory there, which, again, is sort of structured in a similar way where you've got software, fastqc; you have a main script here which is essentially just a workflow that is calling this main script in order for it to be tested and then you have this test YAML which is just a YAML5 file containing, for example, md5sums for the output files generated by fastqc and so, for any given tool these are the files, generally, that you would need to change.

[19:05](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1145)
There's one or two more if you were to submit a pull request to nf-core/modules but it's literally as simple as that and once these have been added (I think that’s the tricky bit) once you add this module to nf-core/modules, it's there and that in the worst-case scenario we may have to change a few md5 sums because things have been updated across releases of it all; but once that tool is there, then we can work with that and that's why I think it's really important to, sort of, fill this out.

[19:33](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1173)
So, this is what a typical module file looks like in our current syntax. This is, obviously, likely to change depending on what gets added to Nextflow and other features that you know or updates we decide to add, but for now this is a typical process where you have just a process name and some sort of publishing logic, containers, your inputs, outputs and a script section. So, I'll try and break this down so apologies for the dense text, it's the only way I could think of to sort of summarize this information to you and take you through it one by one. And these links work by the way, so when I make this presentation available you should be able to just click on these and and get to where you need to be in terms of where this code is on github

So, one of the one of the more important things with (goodness me, 20 minutes) with options, with modules is that we need to be flexible to be able to pass options around two modules and so this is important, for example, you may not always want to publish the file in the same place or you may want to pass different arguments to a command line tool that you're that you're using as a different developer (so you might have installed samtool sort in your pipeline and you want to give it different command line arguments and I may want to use other ones).

[20:58](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1258)
And so, I came up with a simple set of options that you can actually use to do this and, again, these six options that that I've listed here were enough to deal with an end-to-end rnaseq pipeline so they're pretty flexible. It may not be perfect but it works and the idea is that these options are initially initialized here within the module file and then, as I'll show you in the next slide, these options can be overwritten by the parent workflow using the include statement. So, then these options can also be provided to this module file and overwritten from the parent workflow.

[21:35](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1295)
By default, they're initialized to empty values like strings or false, but they can be overwritten and that's the key thing here. You can use those options to then overwrite where you're saving files, for example. You can also provide any non-mandatory arguments as strings to this module. And I think that's quite important in terms of flexibility, so all mandatory arguments (what we would consider mandatory arguments) that should be defined in the module are anything that that involves using inputs and output files because those (in sort of the ethos of Nextflow) need to be staged properly and put in the right place and they need to be defined as inputs and outputs.

[22:14](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1334)
Also, for example, anything where you can use Nextflow resource variables and define them in the script, like threads and so on. So, to show you how these are how these options are passed around (I mean, it took me a while to figure out how best to do this) but every pipeline will have a conf/modules.config and this is specifically within the pipeline repository and this will contain a list of modules that you have, along with arguments and and custom publishing options that that you that you may require.

And in this case, I've just used a simple example where I've got fastqc, I've set the arguments to quiet and this is a non-mandatory argument, it's just a string that can be passed to the module.

[22:56](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1376)
I said publish it here and I've said, also in terms of the files that I want to publish, publish anything with a html dot extension in the top-level directory of my_fastqc and anything with a zip extension, publish it in a zip directory.

And so, it's quite simplistic in the way that it works. Now this modules.config is then typically loaded in your Nextflow config and then all of these parameters then become available to the main script here as a result of this loading. But the great thing about this is now also that because this is a Nextflow params and this modules is just a groovy map, users can overwrite these modules if they want via their own custom configs. I think this is one of the key features here and one thing that I wanted to implement to make things as flexible as possible because with DSL1 pipelines you typically have to physically add a parameter to the pipeline if you want to amend the command line argument, for example.

[23:54](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1434)
With DSL2, you can overwrite arguments now. I'm not saying it's always recommended. Hopefully, nf-core pipelines come with good defaults but in some instances using small or large genomes, indexing may break or some other things may break and so in that case it becomes very useful to pass additional arguments. For example, here I've just appended kmers 10 to fastqc. I've also changed the output directory, again, something that you may not want to do because the pipeline takes care of that but it's possible.

[24:22](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1462)
And also, which files you want to publish so you don't have to publish all of them, you can publish a selected few of them. And then these parameters eventually uh this modules uh instance eventually gets propagated to the main script and then here, this is the key bit here, when you when you use this addParams directory, if you're overwriting the options that I showed you on the previous slide with the ones that you've provided here from fastqc so depending on the combination of what these two configs are you will then provide these options here to the module itself.

[25:02](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1502)
And a real case example is here, where star (which has zillions of parameters possible so I've truncated for good measure there) but you can have all of these parameters provided there.

You can provide one or more of these known types of standardized variables depending on you and the reason for that is because these variables are initialized at the module level and so you won't get errors if they're not initialized via your config.

[25:31](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1531)
In this case, we're setting sensible defaults at the module level which means that you don't have to provide a value for all of them. Another example is because it's just a groovy map you can, for example that star_align process, if you have a parameter in the pipeline that specifically needs to be evaluated for you to add another option to that argument you can just do it because it's appending to a string. If you want to publish files, say if you've got a parameter that says save_unaligned files, then you can just put the additional files within this groovy map and then it will publish those as well.

[26:07](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1567)
It's quite flexible in the way that you can create and pass these options to the modules and as I mentioned before, you just pass them using this addParams directive.

We have had to write a few functions to customize some behavior and also, again, for simplicity in terms of the way that we're dealing with these files. So, there's only three functions that I've had to write that we are now importing into the module script and these are just custom groovy functions mainly used for publishing files and also passing arguments to these command lines tools. So, in this case, as I mentioned before, you are initializing these options here at the module level if they're not provided by the workflow and so they set sensible defaults.

[26:54](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1614)
For example, getSoftwareName here will just use the task.process Nextflow variable and it will just get this bit here. So, by default, this module will publish the files in your output directory and a folder called samtools; which is why in the previous slide I didn't have to provide an output directory for fastqc because I just want the reports to go in a directory called fastqc anyway so i don't need to overwrite that.

[27:20](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1640)
Similarly, these, these and these are all again initiated here in this initOptions map which just initiates decent values now. This works, it worked for the rna-seq pipeline, it may not be as comprehensive, it's slightly annoying where some tools require three arguments and so you have options.args3 which is slightly ugly but it works until we have a better solution.

[27:44](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1664)
And this, here, is how we're saving files again, a couple of Groovy functions and that allows us to provide the files, as I showed you before, in this sort of format where you just have a Groovy map of which files you want to publish.

Not everyone will want to publish all of the files. Some files they may want to publish in different ways.

[28:08](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1688)
Right, so that's just going back to the include statement and where these functions are included in the module. The name of the process is all lowercase, so this is the path to the process, it's samtools/sort. The process name must be the same as the module name, but all uppercase and, again, it should just be separated by a single underscore. This is quite important for standardization, as I mentioned to you before, we're getting the software name, we get the process name, so it'd be nice to sort of standardize this so we always get the main tool name that we can publish as a default.

[28:44](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1724)
And so, this is sort of a niggly thing, I mean there are obviously edge cases where you can have three layers of tools and so on, but this should work for 99.99% of cases.

The tag, by default, uses this meta map of sample information and this is generally provided in the input section here (which I will come back to later) and so this just allows you, when Nextflow is running nicely in the terminal, to see which sample is running and it just gives it a tag as to as to what's running.

[29:14](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1754)
This may not be always possible because in some instances, where you're indexing genomes for example, you don't need any sample information so you just tag it with whatever is appropriate. Appropriately resource labels we have in the conf/base.config of nf-core pipelines and we just have some simple process labels - process low, medium, high, maybe retry or something like that. This just allows us again to just reuse these labels across nf-core pipelines and even in your own Nextflow pipelines if you so wish to use them.

And that just allows us to standardize this a little more in terms of how we're using it. If you don't have that label then nothing will happen anyway.

[30:00](https://youtu.be/ggGGhTMgyHI?t=1800)
The saving files, this actually took me quite a while figuring out how to do this properly but I think we've got something quite simple now. So if you've seen nf-core pipelines, you don't need all of that extensive if-logic and so on it's literally just simplified to this one line which is calling that function `saveFiles` in `functions.nf`. The crux of it is that the `publish_dir`, that you provide via `$options` is just the published directory above the main output directory so it would be `outdir/fastqc`. You can also choose to publish by id, so in this case it would be sample id, for example, and in that way you would have the output generated per sample. You also have this `publish_files`, which I showed you, which is essentially just a map of the file extensions that you want to publish in your workflow and whether you want them to go in a separate directory. This has tended to work quite well; if you don't provide `publish_files`, everything will be published and if you set it to false then nothing will be published.

[31:03](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1863)
One tricky thing that took me a while figuring out is how to actually get this working with Nextflow ´-resume’ because Nextflow only caches the process level stuff, so that the script and and all of that other things it doesn't really cache the the saving and the publishing functionality which is how it should be doing things.

In some instances right early on in the development I found that saving or changing things with saving would break the caching ability but that's just because I was doing things in the wrong way. Eventually I had a light bulb moment, had a conversation with Paolo, he was like, use that params, done that and everything seems to have fixed itself.

[31:41](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1901)
Containers we're using from Biocontainers by default this is why these modules are so self-contained. We're providing these definitions within the module files. We've just realized recently we can't use build ids with Conda and so we can just use tool and version and not, for example, this build id, we can't have that here and the reason is because on different platforms you may have different build ids and that won't work, that will break this module on different build ids and that's something that Anders pointed out a week or two ago.

But as I mentioned, for anything else that is on Conda, we're using Biocontainers directly. These have a build id and they also are mirroring exactly the same singularity image and so you can use these directly within the module file and it’s been working. There's been a few teething issues with installing and using them but overall, it's worked really well.

[32:34](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1954)
So, I'd like to thank Bioconda, Biocontainers for making these available. You can also now build multi-package containers, so if you want samtools and BWA in the same container and so on there are ways to do that. We've used a bunch of those already in pipelines as well.

[32:52](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=1972)
The input (so this is the meta map bit that I said I would come back to), the key concept behind this meta map is that it contains all of the sample information that you need to propagate your channels through that workflow so, in this case, you would have, for example, an id which is your sample name, whether that sample is single end or not (again this is quite genomics focused but I think you get the gist), strandedness of rnaseq (whether it's reverse, forward or unstranded) and we have set some standard ids that we're recognizing (or random variables that we're recognizing) within these module files and so if you have meta id set, for example as I mentioned here in the tag (it's grayed out now) it will recognize that as the sample name.

So, every time this process is run in the terminal it will print the sample name and also whether the sample is single end or not and now this offers another layer flexibility where the module doesn't have to be or the pipeline itself doesn't have to deal with the single endedness or the paired endedness of the sample, the module is now dealing with that in itself so you can have a mixture of single end and paired end samples being provided to the main workflow and the module will take care of that.

[34:05](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=2045)
That was the whole philosophy behind that. Meta, again, may not be required in all instances because you may, for example, just be indexing a genome and any parameters that must be evaluated in the context of the sample must be provided within the process so single and paired and, again, to allow that flexibility we need to add those if else statements for whether things are single end or paired within the module file itself.

[34:27](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=2067)
Hopefully, to deal quite efficiently with file formats and correction and so on, compressed is good. The output if you must emit a name channel. This is useful because for someone else that wants to use your module that you've submitted to nf-core/modules, it helps for them to be able to access the elements that are being produced by this output individually.

[34:54](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=2094)
You might have 10 or 20 different files being produced by the module and so it helps to have a standard named convention for this.

So, if meta is provided as an input, it must be provided as an output and, again, compressed files are great. The script section, as I mentioned earlier, the software name is automatically obtained by this `task.process` Nextflow variable. A command must be provided to get the software version when you submit this module. This has become really useful when we're collating the software versions at the end of a pipeline.

If this is sort of dealt out and propagated to the module file itself, then we don't have to worry about that. It's done once and you forget about it.

[35:39](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=2139)
So, only define input and output files as command line parameters, as I mentioned before, so you you don't define optional arguments. These must come in via this optional $options.args. Similarly, anything that requires (or can use) the Nextflow $task.cpus memory or time has to be, actually time maybe not, but cpus and memory, hopefully in cases where you actually really have to do it because memory gets tricky.

But cpus definitely have to be defined if the tool supports multithreading, for example, and then you have the option to customize the name (output names), we're using this prefix logic here, so you can call your bam file whatever you want if that options.suffix argument instead said .sorted or .mark_duplicates, that will get propagated from your original modules.config to the module here.

[36:36](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=2196)
And so, there's a bunch of stuff to do (sorry, I've been babbling on for quite a while). It's a lot to cover and I was scared this is going to go over, but there's a lot to do, we've added some functionality to nf-core tools for listing, for installing these modules. As I mentioned, Kevin is now doing the linting and md5.

[36:55](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=2215)
We'd also like to be able to create these modules on the fly to update them, to remove them, to check that they're intact and also maybe version check and other stuff.

It'd be nice to have some sort of more automation set up to update the module files themselves because I, and a bunch of others, have physically had to update all 50 or a 100 of these module files at a time until we find our feet with stuff.

[37:20](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=2240)
I mean we're pretty much there but it'll be nice to have sort of a syncing type functionality that we have with pipelines. This has been brought up quite a bit, so the fact that we have to have a `functions.nf` file shipped with every single module file and the fact that that has to be duplicated everywhere; this mainly is that for the module itself to be self-contained. So we can have different versions of the module, the `functions.nf` file itself within the same repository and it won't break things, but we are definitely trying to find a better solution to that.

[37:51](https://youtu.be/ggGGhTMgyHI?t=2270)
More CI - we love CI: reporting linting tests on PR's for reviewers to make it easier; standardize the test data which we're doing already; differences in md5 sums just to help with debugging.

Most, in fact all, of the modules I think at the moment are using biocontainers. It'd be nice to have some that aren't and maybe figure out a way to deal with that.
There’s the workflow package manager, so Junjun Zhang in Toronto, as part of the ICGC Argo are doing something similar. They started out at a similar time where we're just having discussions together about how to make it work and they have sort of branched out and done their own thing with that so we have had brief discussions about how to bring this together but time’s been nuts recently. It’d be nice to revisit that and see where we can compare notes on that.

[38:39](https://www.youtube.com/watch?v=ggGGhTMgyHI&t=2319)
And obviously core Nextflow, as I mentioned before, you can use Groovy and other things to get around a few things but there are a few things that are on a parallel radar, for example dealing with optional inputs.
Optional outputs are completely fine and they work beautifully but optional inputs and where tools can take one or more different inputs, publishing files and so on as I mentioned I've added that customisation to deal with that and it works beautifully but once these are added to core Nextflow we can hopefully strip out some of that stuff there.

[39:11](https://youtu.be/ggGGhTMgyHI?t=2351)
So please come and find us on slack #modules, there's the nf-core modules repository, Twitter, with amazing videos that will be going on Youtube so, again, that's now a persistent resource of information and there's our publication at the bottom.

Thank you to the nf-core community, and the Nextflow communities, Edmund Miller, José, Kevin Menden, Maxime and a bunch of other people that have really helped drive this recently. Thank you all, Paolo as well for his little tips and knowledge and also the Bioconda community for providing all these amazing containers.

And my group for being amazing through this pandemic: we're still standing, getting through work. Thank you to them.

The hackathon is very soon, sign up now for free stash, it's online only: - don't be shy! Thank you!

</details>
