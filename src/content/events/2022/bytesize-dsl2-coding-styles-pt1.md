---
title: 'Bytesize: DSL2 Coding style recommendations (Part 1)'
subtitle: Maxime Garcia - Barncancerfonden, Stockholm, Sweden
type: talk
start_date: '2022-07-05'
start_time: '13:00 CEST'
end_date: '2022-07-05'
end_time: '13:30 CEST'
youtube_embed: https://www.youtube.com/watch?v=KnYPzZ0Dd-Y
location_url:
  - https://www.youtube.com/watch?v=KnYPzZ0Dd-Y
  - https://doi.org/10.6084/m9.figshare.20238969.v1
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: DSL2 Coding style recommendations (Part 1)

This week, Maxime Garcia ([@MaxUlysse](https://github.com/MaxUlysse)) will explain us what he thinks the perfect DSL2 coding style should be.

Come and join the ensuing discussion, it's going to be lively 😉

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=1)
(host) Welcome, everyone. I'm Franziska Bonath, I'm the host today and here is Maxime and he is going to talk about coding styles with DSL2. And off to you.

[0:13](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=13)
I'll share my screen. Hello, everyone, Maxime here. Today I'm doing another talk about DSL2. This time I will try to focus on coding style recommendation. I will not talk much about syntax, but more about organization of the code.

[0:40](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=40)
Just a quick overview. First, what has changed with DSL2? You already know the answer to that, since it's my second topic. What are modules? And so then last, what I think we should do with that. To begin with everything, as usual, this is a disclaimer, these are my own recommendations. I think some other people are agreeing with me on some of these views, but other developers might have other views on that. And we are still trying to forge the best practice and trying to figure out what is easier to read, what is easier to understand. It might - and it probably will - still evolve when we are getting there and I might probably change view on some of these topics. But at the moment, this is what I think we should do.

[1:39](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=99)
What has changed with DSL2? If we follow the Paolo announcement from two years ago on the Nextflow blog, I linked the whole blog post. But for me, the most important phrase is this one: "The module file is nothing more than a Nextflow script containing one or more process definitions that can be imported from another Nextflow script.". This is what has changed with the DSL2 module. As you guessed, I will be talking about modules a lot in this talk. But what actually does that mean? What is a module? So of course, obviously a module is a module. If we follow this definition that would be just so. A process can be a module as well. A subworkflow can be a module and the workflow can be a module. And they all can be interlinked together, which can be a bit confusing, but actually it's pretty clear. But to get even clearer, we agreed at nf-core to have some proper definition. For us at nf-core, a module is a single atomic process that can be called into another script. A subworkflow will be a few chained modules. And then a workflow is an end-to-end pipeline. With that, it's fairly simple. And with this definition, we can decide how we want to organize the code and how we want to do things properly.

[3:21](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=201)
I will explain what I think should be done. For me, the code should be easy to read. That's why it's easy to understand. It's easy to share, easy to modify, and easy to contribute. Because that's what we want. We are all working in open source science. And that's what we want, that's what this community is all about. It's about sharing our code, sharing our work, and working together to achieve these goals. For me, that's what is the most important part.

[3:53](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=233)
Now we are going for some examples. I followed the same bit of code from the module level to the subworkflow level to the workflow level. Just to explain how all should work. The first statement, as I said, all the code for the process is in "modules". That's the nf-core repository. Either we can have a local module within the pipeline or nf-core modules within the nf-core repository. In this case, I'm going to showcase the "ensemblvep" module. The code is fairly simple for a module. We have as usual the tag definitions, the labels that specify the resources that we can decide on afterwards. Then we have the virtual environment or the containers that are specified for that.

[5:02](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=302)
Here we have the input. As usual for the input, first we have the actual file that are depending on the sample that we want to analyze. And then we have the reference file or the reference values that are needed for this, the mandatory one. Then we have all of the optional files. We can specify some outputs. We have the version which is what we want in all of our nf-core modules. Because that way we want to be sure that we can have the possibilities that we want. And then in this case, we have some optional output. Otherwise, you can have some regular output indicator as well. We are doing also a `when` statement in the nf-core module. And then we have the proper part of the script, which is actually just called the tool, and some extra specification to specify some extra arguments or other part from the code. At the end, we just specify the version. This is just the part of the code which is modular because that's what we wanted to do in nf-core. That way we can share the code with everyone.

[6:28](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=388)
As a companion to that, we have some modular setup in the config file. For that particular matter, closures are definitely our best friend. This is what I said in my last talk and this view hasn't changed. With the closure, you can dynamically specify what you want inside. At nf-core, we decided to use custom namespace `ext` directive that allow us to have some specific namespace that we are using. We are using them for the arguments, "args", "args2", "args3". We are using them for the prefix. We can even be crazy and use that for `when`. And then we can also use closure and use other directive at the process level, such as the published year, which is fairly common to use, I guess. If you're feeling very, very crazy, you can even go and change the container.

[7:32](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=452)
This is, for example, what we are doing at the moment in Sarek, still for this module. Within Sarek I currently have a condition to specify if I want this selector to be available or not, because otherwise we have some warnings that this process does not exist. And I just don't like to have a lot of empty warnings. This is optional and I'm hoping we will get rid of that at some point. But this is just starting here. We have a prefix for this specific tool that we want to use. This is what will happen for this prefix. Here we have some arguments. In this case, it's a bit complicated, but basically what we do, we have some basic arguments that will be used in all cases. We also have some specific arguments that will be dependent on the input parameters that we specify on the command line.

[8:44](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=525)
In the end, it's actually fairly simple. We just join all of the arguments together. This is a list. We join all that and then we trim to get just one single string in the end that we can directly put into the process with the `args` directive. Then, because in Sarek we like to do things differently, we specify a specific container in that case. And then we are using our `publishDir` directive to specify where we want to save our file. Depending on the extension of the file, we might have some specific location for that. And that's all. In the end, it might look complicated, but it's actually fairly simple.

[9:38](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=578)
At the subworkflow level, we can actually have several layers of subworkflows. For me, for the first layer of subworkflows, I try to keep it as simple as possible. What I want to do in this level, I just want to chain modules together and just do some tiny channel manipulation if I have to. For example, you need to remap the output from one module to go into another module. Then yes, you can do that at this level. It will, for example, arrive to that. This is the subworkflow that I'm using, that we are using in Sarek to call the module "ensemblvep" and then call the "tabix" module to `tabix index` the .vcf file. This is what we are doing. As usual, we begin the workflow by taking the input data that is related to the sample and then all of the reference genome and optional values that we need to share.

[10:51](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=651)
Then here, fairly simple, still I call the first module. I call the second module on the output from the first one, and I `emit` everything out. And of course, we gather all of the versions for the tools that we use. That's still at this case. If I want to go still, we are at subworkflow, but then it's a fairly higher level because the subworkflow, we can still include other subworkflows in the subworkflow. What we can do there is that we can chain modules together, or we can even chain subworkflows together, or subworkflows and modules. We can also manipulate channels. What we can do here but is not good to do at the previous level, is to specify some execution logic with an `if` block. It will look like that. In this case, I'm calling three different subworkflows. Actually, I'm just calling two different subworkflows. I'm calling the "ensemblvep" that I just showed, the "snpeff" subworkflow. And I'm calling the "ensemblvep" twice because one time I want to use it as it is, and one time I want to use it on the output of the other subworkflow. I need to do an alias to actually be able to use it twice in the same subworkflow. As usual, it takes as an input the files that are related to the sample, and then the reference data and the other values. Then, as I explained earlier, I'm having this `if` statement to control the execution of the subworkflow, and then I collect all the files that need to be collected. Same thing for that. Same thing for that. And then we `emit` everything back.

[12:56](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=776)
If we follow this logic, we can get some fairly easy to read and easy to understand organization of the modules and subworkflows. That way, it's easy to understand where to contribute, what to do, how to change, and how to evolve stuff.

[13:18](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=798)
At the workflow level, it's where we should do everything else. We can still, at the workflow level, call a single module. For example, you might want to call just the MultiQC module here at the workflow level. Of course, at the workflow level, we want to chain several subworkflows, because if we don't do that here, where are we going to do that? And then, of course, you still want to do some channel manipulation, because yes, that's your main workflow script. That's where you want to do all of the magic. And, of course, the execution logic still happens over there.

[13:57](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=837)
This is what happens within Sarek at the moment. Still here, I have my execution logic. If my input parameters are right, then I'm going to do that. And here, I just call my subworkflow within the main workflow and I gather all of the used software version in the report that I need to add. That's all. For me, this is how we should organize our code for the module, depending on if we are at the module level, at the subworkflow level, or at the workflow level.

[14:43](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=883)
I also have some small syntax recommendation. We are trying to make proper recommendation guidelines on the DSL2 syntax. We are working on the document with several other people from nf-core. If you want to contribute, the link is in the title here. And what I would to say is, first, indentation is your friend. That's a fairly good statement. And a lot of people that are coming into bioinformatics will learn to code with Python. In these cases indentation is already deeply ingrained into your habit. Let's continue working on that and let's keep indenting. It makes the code nice to look at. And I like to have a nice code to look at.

[15:33](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=933)
In a process, I saw several different ways to collect several files, just to call it with the same parameters. And this is a proper way to do so. I try to enforce that in all of the GATK4 modules, but that might have escaped me at some point. This is fairly small and it's a small one-liner that is easy to understand. Then for the channel assignments, that's something that you want to do in a subworkflow or in a workflow. I personally prefer to specify which channel I want to assign things to first. But some other people might prefer to have the channel in the end. I personally prefer the first version. I can see that some people coming from an R world would prefer the second one. We really need to figure out with the community what we want to do. But I'm going towards the first version of this line.

[16:48](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=1008)
Some last tips before I open the discussion. Please add comments to your code. It's good for you. Good for your colleagues or anyone who are going to have a look at the code. It's also good for future-you because, yes, believe me, two weeks from now or three months from now, you're going to look again at your code and you're not going to remember why you did that. You want some comments and right now you're your own best friend. You should be able to help yourself and do that.

[17:21](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=1041)
Also, it might be important, I notice this sometimes within nextflow, it might be difficult to differentiate between queue and value channels. And it could be a good idea to have that difference in mind from the beginning when you design your pipeline. Usually we do use value channels for all the reference files. And sometimes you don't understand why your process is not executed several times. It's because this channel that you believe was a value channel is actually a queue channel. We need to be careful with that. Otherwise, good a tip would be to not hesitate to ask questions. We have Slack. We are available on that. And I personally like to discuss with other people and stuff. Discuss, that's good. If you're on site with other people, have a coffee break. That's always good to start discussion with anyone else.

[18:28](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=1108)
I'd to thank all of the institutions I work for and all of the institutions I work for and with on my project. I'd to thank all of the institutions that are working with us on nf-core. And I would to thank all of the contributors that we had so far at the moment on nf-core as well. If you need some help, of course, Slack as usual or everything else. If you need more help, we have some documentation and we have some previous bytesize that you can check up on. Otherwise, I'm open for questions now.

[19:09](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=1149)
(host) Thank you very much... I don't seem to show... Anyway, I have now opened the possibility for everyone to unmute themselves, if you have any questions. Otherwise, we have also questions in the chat.
(speaker) OK, let's begin with a question in the chat.

[19:30](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=1180)
(question) I do have a question.
(host) Oh, yeah, go ahead.
(question cont.) OK, so I guess this is my first question, which is I'm wondering which specific Slack channel discusses Sarek in the nf-core workspace. Because I haven't found any specific one for Sarek or it's just one of the existing channels.
(answer) No, no, we have a specific Sarek channel to discuss about Sarek.
(question cont.) OK.
(answer cont) We have a specific channel for each pipeline and we do have a specific channel for every main topic. Otherwise, whenever you don't know anything, don't hesitate to go to the #help channel and then we will direct you towards the right channel.
(question cont) OK, so it's #sarek or something that, right?
(answer cont) Yes, #sarek.
(question cont) Oh, I see. OK. Oh, I found you. Thank you.

[20:32](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=1232)
(question) My other one is more of a comment as a follow-up to what somebody said about the last comment about a user of a pipeline. Yes, so DSL1 or DSL2 are not really important. However, as a developer, it matters. Yes, so DSL1 to DSL2 conversion is a lot of work. I tend to disagree with the first part of that, which is for the user, DSL1 or DSL2 is not really important. It _is_ important because people are forced to actually go back to their code for those that actually did DSL1 without specifying how their code should run. Since nextflow pretty much just force people by defaulting to DSL2 and the code breaks. People are now forced to revisit their work. If it wasn't that important, even for a user, then that would not be necessary. I would tend to disagree and say it is, in fact, important to have pipelines in DSL2, whether it works in DSL1 or not is irrelevant at this point. Because unless if, you know, you were thinking about the fact that DSL2 would be default when you were coding and then you had made it very clear in each of your scripts that it was DSL1 enabled and stuff like that. Right. I just wanted to point that out. But thanks. That was a great talk.
(speaker) Thank you.

[21:57](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=1317)
(host) Yes. Also a question from Matthias.
(speaker) Okay, let's go.
(question) Yes, thank you for your great talk. Could you please elaborate a bit more on the extra arguments, because I noticed, for example, that you've wrapped them in brackets in your module example.
(answer) Yes, let me go back to that.
(question cont.) I haven't seen that so far.
(speaker) That was that or not at all? Here then? In the module?
(question cont.) No, it was part of the code where you had the...
(speaker) Not this one either, right? Sorry, I'm trying to find it back. No, not this one.
(question cont.) It must be the config where you provide the extra arguments.
(speaker) Are you saying the for the `publishDir` thingy?
(question cont.) Exactly here. Exactly here. The param stopped that log tree or something. You've wrapped them here in brackets.
(speaker) Oh, yes. Here what we're doing, so we have everything all of this `ext.arg` is in brackets. That's a whole list.
(question cont.) Yeah, that's the list. And then you... Not brackets, braces.
(speaker) No worries, it's fine. I always confuse the name of this, even in French. Here it's a list and for each element of this list, I specify either directly the params or directly here I have a tertiary `if` (which is a bit long, sorry), which let me specify, if (sorry, I'm going back there)... If this statement is true, then I have this first string that will be considered. Otherwise, it will be the second string, which is an empty string.
(question cont.) And you basically just have the braces for consistency. If you don't evaluate the end, then it doesn't matter whether it's in braces or not. Or is that...
(speaker) No, no, the brackets are important because otherwise if we don't have that, then it's not a list and I want to have several arguments. I could make one complicated if, getting all of the different stuff, but what I want one string on which I can append other arguments on top of that.
(question cont.) This is completely clear. I'm just wondering about the limitations because sometimes it's you've wrapped curly braces still around. This is when you have a meta map in there.
(speaker) Ah, okay, you mean that here?
(question cont.) No, probably we should just take this.
(speaker) Yes, let's take this somewhere else. There should be some easy stuff for that.
(answer) He's talking about when to use a closure. The closure for example is used in `ext.prefix`. When you use a closure, the variables are evaluated during the task execution. But if you don't, so for example, the `ext.args`, this one doesn't have the curly braces, so this is not within a closure, and you can't use any variables from the task execution context. And it's evaluated as soon as the config is loaded. Right at the beginning of the workflow before any of the tasks are executed, but if you use a closure, then, one, you can access variables within the task context, anything from input. But also, it means that things params out there and their evaluation is delayed until the execution of the task.
(speaker) Thank you Mahesh. No, I understand what was the question. Sorry, I misunderstand you.
(question cont.) I didn't make myself clear. But it's still a big mystery to me in this regard here. The details.
(speaker) Then don't worry, we can have a more detailed talk about that another day.

[27:00](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=1620)
(host) There's a question from Phil.
(comment) Hi everyone, thanks for the talk Maxime. I just wanted to reiterate something that came up in the comments on zoom just now, mostly in case anyone's watching this on YouTube at a later date. It was a really good comment about DSL1 and 2, and how running a DSL1 pipeline with newer versions of nextflow will crash in a slightly nasty way. But just to note that that's actually only a fairly specific version of nextflow where that happens. Paolo switched it to DSL2 by default and we saw these crashes happening. We spoke to him and now the newer versions of nextflow since version 22.04.3 should automatically detect whether the workflow is DSL1 or DSL2, so you can go back to just running it without any flags and it should just work whether it's DSL1 or whether it's DSL2, no nasty crashes anymore.

[27:57](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=1677)
(question) But will we keep DSL1 inside the main nextflow or at some point will we completely discontinue it?
(answer) At some point it will be completely discontinued. In the future... it's some time... it's mentioned in the nextflow blog. I think it's planned for kind of mid 2023, the released versions of nextflow will stop supporting DSL1. And from that point, any DSL1 pipelines will have to be run with older versions of nextflow.
(answer cont.) Yes, but then it's still fine because we can still install older version.
(answer cont.) Exactly, so the workflows will still run just not with latest version of nextflow. And yeah, if you're interested in converting workflows into DSL2, then we have a slack channel, an nf-core dedicated to this topic called DSL2-transition. And the same goes also for subworkflows.
(answer cont.) Thank you, Phil.

[28:59](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=1739)
(host) There's another question from Olaitan.
(question) Okay, thanks. Okay, good. I'm going to make it brief. My actual question is about Sarek again, the reference genome has to be one of the existing ones, maybe GRCh38, 37 and the likes, right? So if I have a scientist, who has his own reference genome, which is not one of the standard types, how do I handle the situation?
(answer) Then it's fairly simple. We have these parameters that you can choose within Sarek. But can we talk about that more in the Sarek channel?
(question cont) For sure, this is going to be a long talk, for sure. Thanks.
(answer cont) No worries. Just don't hesitate to ping me on slack and I will answer directly.
(question cont) Perfect. Thank you.
(speaker) Thank you.

[29:53](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=1793)
(host) Okay, are there any more questions from the audience? There's a question for... Ah, no. Phil also posted a link to Nextflow.io?
(comment) It's the Nextflow blog post talking about the end of DSL1 support. I'll put it in slack as well.
(host) Thank you. Also earlier Mahesh posted a link to Carpentries training material for coding practices. I will also leave this link in the slack channel for bytesize.
(speaker) Maybe I can include that in my slides as well.

[30:26](https://www.youtube.com/watch?v=KnYPzZ0Dd-Y&t=1826)
(host) Yes. Okay, if there are no more questions from the audience. Thank you so much Maxime for the talk. I also would to thank the Chan Zuckerberg Initiative for funding these talks. And as always, you can continue these talks in slack on those hashtag bytesize, or specific to any of the channels for different workflows, and thank you very much.

</details>
