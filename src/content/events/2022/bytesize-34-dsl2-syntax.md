---
title: 'Bytesize 34: Updates on the new DSL2 syntax'
subtitle: Maxime Garcia - SciLifeLab / Karolinska Institutet, Sweden
type: talk
start_date: '2022-02-22'
start_time: '13:00 CET'
end_date: '2022-02-22'
end_time: '13:30 CET'
youtube_embed: https://youtu.be/17NqUsh73BU
location_url:
  - https://youtu.be/17NqUsh73BU
  - https://doi.org/10.6084/m9.figshare.19263401.v1
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 34: Updates on the new DSL2 syntax

This week, Maxime Garcia ([@maxulysse](https://github.com/maxulysse)) will tell us all about the latest changes to the DSL2 template.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://youtu.be/17NqUsh73BU&t=1)
(host) Hello everyone, thank you for joining us for yet another bytesize talk on this Tuesday. I would to begin by thanking our funders, the Chan Zuckerberg Initiative for supporting all nf-core events. A few details, the talk will be recorded and is currently being recorded and the video will be shared on our YouTube platform and shared on Slack as well. For today we'll have as usual a 15 minute talk and then it will be followed by a question and answer section where you are free to post your questions on the chat box or unmute yourself and ask your question to the speaker. Today we are glad to have with us Maxime Garcia who is a bioinformatician at the Science for Life Laboratory at the Karoliska Institute in Sweden and he'll be discussing or talking to us about DSL2 which is a syntax extension that implements the definition of module libraries and also a better way of writing more complex data analysis pipelines using Nextflow. Maxime, over to you.

[1:07](https://youtu.be/17NqUsh73BU&t=67)
Thanks a lot Simeon. Hi everyone, Maxime here. Today I'm going to talk to you about the new update that we have on the new DSL2 syntax especially for the modules. Brief overview about what I'm going to talk about, what's new, what can be done and what should be done. First let's begin with a disclaimer. This is my own takes on the new syntax. Other developers might have some other idea and the best thing about nf-core is that we are a community. Of course there are some driving forces such as Harshil or Mahesh that try to do stuff and that are having a lot of ideas and stuff. But what I like is that everyone has their own voice and everyone might say and might improve everything, even me. Of course because we are all doing this as a community the syntax and/or the logic will definitely evolve. This is just the current state of the new DSL2 syntax. What we want to do is, we want to forge the best practices for all that. I will just show you what I'm doing right now and let's discuss at the end to see where this is going.

[2:36](https://youtu.be/17NqUsh73BU&t=156)
What is new in the new DSL2 syntax for the module? Now the modules are fully self-contained. We don't need any function.nf anymore and we don't need the params whenever we call a module sub-workflow. All of the logic when calling a module or sub-workflow can be done with using the new task.ext directives. We can set up a different argument. We can set up the `prefix` for the file name and we can use the `when` directive to set up when and if a particular module should be run or not. Of course this task.ext directives are now used with a `withName` selector in the module.config. Instead of using a huge params.map, now we just have a `withName` selector. It's a brand new world. We can do almost everything and that's the main issue with that. Because we can literally do almost everything and with that the thing that we really need to figure out is, to make good closures, to decide how to use the argument when we have several or just one, we can use closure to decide that. We can use closure also to design the prefix and we can use closure to decide when to run the module.

[4:17](https://youtu.be/17NqUsh73BU&t=257)
But there are some downside. We've got this new syntax for the module that allows us to use the task.ext directives. That means that the logic can go into the config. The main issue is that now we have the logic in the workflow, in the subworkflow and also in the config, which can be messy. That can be super bad. What should be done is that we should be careful whenever we are setting up task.args, especially because for the argument, if we are setting that up in the config, then we must be very certain to explain how it's happening and where and why and everything. I advise to everyone to write comments to explain the whole logic of that. First it will be good for future you because you will forget about it. I know I already did some stuff so the comments are already helping me and I just started last week. Definitely it's good for everyone and it's good also for other developers because as I said earlier we are part of a community and we are not coding only for ourselves. We are coding for everyone. It's good for everyone and comments are always a good practice.

[5:39](https://youtu.be/17NqUsh73BU&t=339)
Now it's time to have some examples. I will begin examples with Sarek because I just finished a PR that will merge with Friederike last week and we are now the leading pipeline with all of this recent development. If I look at the prepare_genome subworkflow which will prepare the indices that we need and will prepare some other files that we need before launching the whole pipeline. It's now super simple. I just launch all the tools and that's all. The logic behind all that will be in the config file. As we can see the subworkflow is super simple. It looks super clear. It looks amazing and I'm very happy with that. Which is why we need to be careful and we need to comment. Here for example I commented. I said that this will be run if the aligner is bwa-mem and this one is run if the aligner is bwa-mem2. Of course I added some specific information here in the beginning of the file, to set that for all modules, when the close condition is defined in the module.config, to determine if the module should be run. Here I explained how the condition is defined in an extra comment. To say if there is an extra condition then it is specified in the comment which is what I just explained here.

[7:24](https://youtu.be/17NqUsh73BU&t=444)
Now let's just have a look at the module.config regarding that. For example if here with bwa-mem1. This is the published directive that will figure out how to publish and save the file, if you want to save the file or not and we can do all that within the configuration. We don't have to take care of that anywhere else and that's so simple. For example here in this case we will save this file only if we specify the `save_reference` params using the `publishDir` mode and with this specific path and the specific pattern. We will run this process only if we have the `param.aligner` which is bwa-mem. Only if we don't have the bwa params. That means that we don't have any bwa indices that are provided to the pipeline and only if we start the pipeline with the mapping step. If we start the pipeline with a later step for example if we start the pipeline with variant calling then we don't need to have the bwa indices so we don't need to run this one. Similar for the bwa-mem2 process we only run it if we have the bwa-mem2 and it's so on for all of the other processes here in this pipeline.

[9:01](https://youtu.be/17NqUsh73BU&t=541)
For the indices and the preparation of the indices and all the other tools it's fairly simple. I just choose some condition within this closure to decide why it should be run and how and why and if and not.

(host) Maxime sorry to interrupt. Somebody asked if you can increase the font size a little bit.

(speaker) Oh yes of course sorry. Like that should be better then.

(host) Yeah looks good to me.

(speaker) Yes, okay. Then let's see something a tiny bit more complicated. Here we will be looking at the mappings of workflows that we use in Sarek and that I hope to publish one day in the nf-core repo so that it can be used by other pipelines as well. This workflow has been refactored several times by Friederike and I, and I'm pretty sure we have other people that are looking into that as well and that we will improve that again. But I'm always happy for that so I think it's good. Here it's the same, in the Sarek workflow we have the whole logic that decides if this subworkflow is run or not. I will not show that, but I will just show here inside the workflow or inside the subworkflow how it goes. Basically with this task.ext directive we can already set up the whole logic inside the config file and so the pipeline itself is much simpler and here we run just `BWAMEM1_MEM` or `BWAMEM2_MEM` on the input file with the indices and we set up "true", because we want the output file to be sorted. Here we're just gathering the bam file outside and we are remapping and we don't want to start the workflow, but that's an extra step and in the end that's all. Only if this closes true we will merge all the bam files and we will only do that if we want to skip the `markduplicates`. It's all explained here in the command: only if we want to skip `markduplicates` or only if we want to save the bam file. Only in this step we will merge the bam files and we will index them.

[11:45](https://youtu.be/17NqUsh73BU&t=705)
Then of course we gather all the versions: Here we have all the modules that are called and the whole logic will happen again in the workflow, in the config file. Up here we see that similarly to what we've done with the indices, we run this bwa_mem only if we have the `params.aligner` equal `bwa-mem`. For bwa_mem2, it's only if we have the `params.aligner` equal `bwa-mem2`. We can see that we set up a particular argument depending on the meta map. In our case in Sarek we have a specific handicap if we have some tumor samples, so if our status is one (meaning it's a tumor) then we have this particular value. Otherwise it's the regular parameters that we use. Similarly we have a particular prefix that we use only if we split the fastq file at the beginning. As I explained with the merge and the mapping, we only do the mapping when we save the bam file, the mapped bam file, or when we skip `markduplicates`. That's all. This whole idea about improving this whole syntax, really allows us to make the subworkflow easier to read but in the meantime we really have to push everything into the comments to explain all that. I made a bad copy paste here because this is exactly the same...

[13:56](https://youtu.be/17NqUsh73BU&t=836)
Let's see one more complicated subworkflow and I think that will be my last example for today. This will be the `markduplicates` subworkflow which can be skipped as I explained earlier. As an input we take the `bam_mapped` which contain the meta map plus the bam file, or the `bam_indexed` which contain the meta map plus the bam and the index. We only have one of those depending if we are skipping `markduplicates` or not. In our case... Oh we can have both of them because it's an optional channel, but let's... it doesn't really matter here. In this case I will run samtools on the bam file to convert the bam file to cram when we have no duplicates, which is why I have this huge name for the module. This will only be run when we are skipping `markduplicates`. Otherwise when we are running `markduplicates`, if we are running `markduplicates` with SPARK we will run that. If we want to have some quality control tool run out of `markduplicates` then the output of `markduplicates` will be bam, otherwise it's cram and if we have a bam then we will convert the bam to cram because we want to use cram in our pipeline. This part is slightly difficult to understand which is why I try to comment everything and which is why I try to put extra comments in the config as well.

[16:06](https://youtu.be/17NqUsh73BU&t=966)
If we are not running SPARK, but run `markduplicates` then we are running the regular module for that, which is the `GATK4_MARKDUPLICATES`. Then we are converting the bam file to cram and then in the end this channel `cram_markduplicates` will contain only one of the following channels, because we only have one solution: either we are skipping `markduplicates`, or we are running `markduplicates` SPARK with bam output or we are running `markduplicates` SPARK with cram output or running the regular `markduplicates`, which is just one of these solutions. In the end, if we are running `markduplicates` SPARK and the report on the bam file then it runs this one and otherwise we run the report on the `markduplicates` bam output or input. Otherwise we do `SAMTOOLS_STATS` on the cram file and that's all. In the subworkflow it's a bit complicated but it looks clear to read and I think that makes it easier to understand even if the logic is a bit fuzzy, which is why we have everything here in the module. Here similarly we will have the prefix that will explain what the output file should look like and we have proper `when` directives that will explain to us how to run it and why and where. This is all explained there.

[17:55](https://youtu.be/17NqUsh73BU&t=1075)
Here we have that for all of the processes that we have there. What can be done here to improve will be to sort out all of the `withName` selectors and I think it was a good idea to first group the selectors by subworkflow, but maybe sorting out the selectors will be good. I'm still thinking about if we should sort them alphabetically or if we should sort them in the order that they are in the subworkflow. That is a different solution, I'm not sure what to do there. Of course what you can do as well with the `withName` selector, you can group several modules together, which can lead to extra issues because you might not notice that you're defining twice the same `ext.suffix` or `ext.when`. You need to be really careful when you're changing several modules at once and that was all for my examples. I want to thank all of my institutes and all of the institute I'm working with and everyone that helps us with Sarek. All of the institutes that are part of nf-core and all the people that are contributing to nf-core. If you need help I will recommend to watch the old bytesize, even if they are not up to date. Otherwise, you know where to join us on slack or on twitter, on everything, and now I'm open for question. I think I saw there was some raised arms.

[19:47](https://youtu.be/17NqUsh73BU&t=)
(host) I think Moritz has a question.

(question) Thank you for the introduction here. I didn't follow all the discussions on slack and github, can you say why this `when` syntax in the config was chosen? Because what I saw now was that mostly in your configs you were referring to global parameters and then in the subworkflow you had a comment of the condition. To me it is much more obvious to put in the subworkflow an if statement with those global parameters because then the logic is right there to read and it hides away that logic and makes it more difficult in in my opinion. But I haven't read the whole discussion around it.

(answer) I agree with you that it hides away the the logic, definitely. The logic is a bit... I don't know why but I think it's a good way to go because for me it will be much easier to control the stuff, how it works, and it will be much easier to make some subworkflow. That will be easily shareable between different pipelines, which I think is something that would really advance at the nf-core level. For now what we are doing is that we are getting good at having subworkflows that are looking good. For example we have a good a tringle or fastqc subworkflow, but that is mainly just copied over from one pipeline to another and I think it will be very good if we can do that and I agree the logic is hidden. But if we explain it well with comments it will be good and you don't need actually to use the if logic. The if logic or the whole logic with the `when` prefix can be decided on in the config or not. That's something that you can decide for yourself in your own pipeline or in your own subworkflow. For me, adding this `when` directive to the module gives us the possibility to do more stuff. The problem with that is, yes it can be good or bad depending on what we do with it. I hope I reply to your question.

(question cont.) Yes, thank you.

[22:30](https://youtu.be/17NqUsh73BU&t=1350)
(question) Frederick asks whether dividing the configs for subworkflow can reduce size and increase the [inadible]

(answer) Yes. Having just a simple config file for each subworkflow could be easier and we could even have the config file sitting in the same folder as the subworkflow. That's something that we can decide or not, I think. That's what I said earlier, we need to decide what are the best practices and how to reinforce that and what to do. How to follow and how to go on with that.

[23:08](https://youtu.be/17NqUsh73BU&t=1380)
(host) I don't know if anyone has another question. Apparently there's no one else who has a question.

(speaker) Okay, then I'm pretty sure we will have more questions on slack as soon as more people have seen that and as soon as more people realize what we can do with that, because definitely this new syntax can be very helpful - or could be very dangerous depending on what you want to do. Especially with this new usage of the `when` directive.

(host) But maybe, even if you develop a standard syntax for the normal processes modules, quality control and trimming. Then they can literally be applied to everywhere, you probably don't have a problem with it. Thank you guys for joining I'll see you next week for another bytesize.

</details>
