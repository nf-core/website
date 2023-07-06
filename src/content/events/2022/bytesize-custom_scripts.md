---
title: 'Bytesize: using custom scripts in Nextflow pipelines'
subtitle: Chris Hakkaart, Seqera labs
type: talk
start_date: '2022-11-15'
start_time: '13:00 CET'
end_date: '2022-11-15'
end_time: '13:30 CET'
youtube_embed: https://www.youtube.com/watch?v=3aA5-s8PAF0
location_url:
  - https://doi.org/10.6084/m9.figshare.21572238.v1
  - https://www.youtube.com/watch?v=3aA5-s8PAF0
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: using custom scripts in Nextflow pipelines

This week, Chris Hakkaart ([@christopher-hakkaart](https://github.com/christopher-hakkaart)) will show how to integrate custom scripts, such as R or Perl o Python scripts, can be integrated into a Nextflow pipeline.

<details markdown="1"><summary>Video transcription</summary>
**Note: The content has been edited for reader-friendliness**

[0:01](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=1)
(host) Hello, everyone. My name is Franziska Bonath. I'm today's host. And with us is Chris Hakkaart. And he is talking about how to implement custom scripts into your Nextflow pipeline. Off to you.

[0:16](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=16)
Awesome. Thank you. And thank you for the introduction. Today I'll be talking about custom scripts and how you can add them to your Nextflow pipeline. I think this talk was inspired by questions that we see on Slack occasionally where people are having trouble implementing a Python script or an R script or a Perl script or some other script as a part of the pipeline. And quite often it's not necessarily mistakes, but things that can be done to make pipelines a little bit better and more readable. Today I'm just going to try and outline some of these things and hopefully everyone will be able to walk away with a little more understanding of how to do this with Nextflow.

[0:57](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=57)
Today I'll be talking a little bit of background, kind of outlining the problems and how we can solve them. I'll introduce my first pipeline, which is a Nextflow script that I've very quickly written to demonstrate how to use the bin and templates directory to store your custom scripts. I'll quickly talk about managing dependencies and some of the things that you might need to consider when you are packaging those together, as well as a really quick summary.

[1:27](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=87)
Background. I don't think it's a secret that in real world pipelines you often need to use custom scripts, which can be written in different languages, so Bash, R, Python, Perl, as well as others. With Nextflow you can integrate any scripting language into your workflow by adding the corresponding shebang to the code blocks or the script block, and I'll demonstrate this really quickly in the next couple of slides. You can avoid keeping large code blocks in your main workflow by executing them as custom scripts. Some scripts can be really short and others can be really long. In either case, to improve readability it's recommended that you store these elsewhere and then execute them using Nextflow rather than just having a big really troublesome code block, which can be quite difficult to get through if you are trying to read through someone else's code.

[2:20](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=140)
This is my first pipeline, which is a Nextflow pipeline that I've written. As you can see it contains a single process called myScript, which is going to take string values as an input and give standard output. In the script block there you'll see that all it's doing is taking the strings defined with the name "str" from my input, and it's going to turn all the lowercase letters into uppercase letters, so nothing too complicated. This is a really simple, in this case, single line of code. In reality, your script could be much, much larger and also written in a different language. In the workflow block down the bottom here, we have "this", "that", and "other" being the three string inputs that we're taking from the channel. I'm pipeing that into my process and then just using the view operator to show that my output in my screen terminal window.

[3:12](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=192)
As I mentioned earlier, you can just easily add a shebang to the top of the script block to change the language of that script block. In this example here, I've changed the script from bash to an R script. And what I've done here is just rewritten my job or what I was trying to do, which is turning lowercase letters to uppercase letters using this `toupper` function, which is a basefunction in R, and then just catting that. Now I'm printing it out to the screen. Nothing else has changed in this pipeline apart from this here in the script block. I just want to point out that shebang again, because in this case, it's an R script, but you could also include Python, Perl, like I said, any other shebang to decide which type of scripting language you would like to use.

[4:08](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=248)
This is just a really short animation of what this pipeline is actually doing. And as you can see, it's just running the pipeline and printing "that", "other" and "this", all in capital letters. This is just those three strings that I've included there being printed out to the command window. Just once again. So it's just executing that pipeline and printing those out.

[4:31](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=271)
While this is a short, single line of code, if it was much larger, you might automatically think, OK, I want to remove this and have it as an executable R script somewhere else in the pipeline or someone else on your system. And this is what I've done here. I've decided to call my code block `myfirstscript.r` just because it's an R script. I have changed this slightly, which I'm using command arguments to allow for an input to be taken as the script is executed in my script block. That's what's happening up here. We've got my first script taking the command arguments. True, so it's going to take the trailing arguments and absorb them as a part of the script. Down here in my pipeline, what we'll see has changed is that we've got this full path to "myfirstscript". And then it's using this "str", which is the input or the named input as a part of this process.

[5:30](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=330)
So what I've done here, and this will work on my system, assuming that full path to my first script is a real file path. This will work on my system, but it's not overly portable. With Nextflow, I think one of its greatest strengths is the portability of the pipelines. There are other ways to do this, if you were to try and share this pipeline with someone else, you wouldn't need to hard code in this file path again. You would just be able to do this automatically using bin or templates, which I'm about to talk about.

[6:05](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=365)
One thing I did want to point out here is that if you are making a script from your code block like I have done here, you do need to make sure it's executable, so you need to run this on your code to make sure that's got the right permissions.

[6:18](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=378)
The first way that you can store your scripts is using a bin directory. So instead of trying to include my first R script as a part of the directory, the same directory as your pipeline, what you can do is create a folder called "bin" and then store your scripts in there. What Nextflow does is whenever you execute a pipeline, it will look for the bin folder within the directory of your pipeline. And if it's there, it'll mount the files or the scripts in that folder to your path, and they'll be automatically executable in your pipeline.

[6:58](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=418)
What that will look like is here: you've got "myfirstscript" in the bin directory. This is in the same directory as you'd expect your pipeline to be in. This is the same script as before, the only difference here is that you're not having to specify the whole file path. You can just have it here in the bin directory. In your script block, the only thing that's changed is instead of having your whole file path, you've just got `myfirstscript.r`. What Nextflow will do is automatically, like I said, mount this file, mount this script in your bin, and it'll be executable automatically, which is a really powerful way of storing the script. You can have lots of different scripts in a bin directory, and they can all be executable automatically.

[7:41](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=461)
There's also another way to do this, which is using the templates directory. This is very similar to a bin directory in that you can have a folder called templates next to your pipeline file with your script stored in there. Now there are a couple of differences that you might notice. While of course it is in a different folder, so it's in this templates directory, here I don't need to specify arguments or the arguments command that you saw previously, because what Nextflow will do is treat this exactly like you'd expect a script block to be specified. And what it will do is, using this template, it'll just look in the templates folder and execute this as if it was a code block included here in the speechmarks. What you'll also notice is that here I've included this named input, and it will automatically be able to use this straight away. You don't need to use arguments like you did with the bin directory.

[8:45](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=1)
The next thing is dependencies. The scripts I have been executing work because I have R installed on my system locally. But if you are using it on a different system or you want it to be 100% reproducible on a different system, you will need to consider dependencies and how they're managed. Dependencies with a custom script are managed in the same way as other tools, and I'll show this very shortly with some examples. But of course with a custom script like this, you might expect to have one or more different tools or packages, which can add complexity to how these are integrated and stored. As with other tools, if you are using multiple tools in the same module or same process, you can store them in a combined mulled container. While I won't go into it today, there are helper tools and documentation available. These slides will be available, and both of these are clickable links where you can read a little bit more about this, about how nf-core has a module's mulled function that can help you find a container with the dependencies you're looking for, as well as how to package multiple containers in one mulled container.

[10:00](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=600)
This is an example of a mulled container.. I mean of how to use dependencies with a custom R script. This is directly from the RNA-seq pipeline. The process is called `salmon_summarizedexperiment`. What you will notice is that down here in the code block, we have `salmon_summarizedexperiment.r`. This is an executable R script with a couple of argument inputs. And at the top here, we have the conda and container declarations. In this situation, it is just a single R packaging tool. This is R base as a part of this package already. You only need to specify this once. And this also already has the Galaxy project and biocontainers, images, containers available. This is a relatively simple example with just the one tool. But you can also have an absolute monster.

[10:59](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=659)
Also from the RNA-seq pipeline here, we have the `deseq2_qc` process. What you can see here is that we have a very large number of different tools that are implemented as a part of this. In this case, a mulled container has been created, which contains all of these. This is probably not a great example because the versions of these tools haven't been pinned to the conda tools. And that's because there were conflicts as this is getting created. But normally, you might expect to see some version numbers after each of these. Again, this is probably a bit more of a monster of a script. I haven't shown this here, but you can find it by going and looking into the RNA-seq repository. And what you can see here is that this is an executable script again with a number of different inputs that can be taken as arguments.

[11:55](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=715)
Just to summarize, what I've covered today is that NextFlow can use custom scripts written from many different languages. Scripts can be stored in both the bin or the templates directory. And both of these will be available to the Nextflow script. Meaning that you don't need to specify an absolute or a relative path as you're executing a script. It's really fantastic to do this because it makes your scripts much more portable and usable by others. Dependencies can be managed using conda and containers. In the examples I've shown, you can see that it can be quite simple or much more complex with the use of mulled containers to help you store all those together with singularity and conda images. And with that, I will finish. I think we're probably about where I thought we'd be for time. And I'll finish on this. And if there are any questions, I'm happy to do my best to answer them. Thanks very much.

[12:49](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=769)
(host) Thank you very much. So anyone can now unmute themselves. If there are any questions, you can just ask them or put them in the chat and I will read them out. There is one question in the chat.

(comment) If you're on a newer Nextflow version... Oh no, that's I think a comment. If you're on a newer Nextflow version, there's also a link to module binaries.

(speaker) Yes, exactly. So I've included probably some simple examples here using the bin and templates directory. You can also store these scripts another way along with modules. But I haven't gone into that as much. There is documentation on this on the Nextflow website. But today I just wanted to focus on probably what I think are the more easy examples, the simple examples, but you probably don't need to do the more complicated systems with complicated techniques.

[13:51](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=831)
(host) Okay. John is potentially asking a question.

(question) Yes, can you hear?

(host) Yeah, we can hear you.

(question cont.) Great talk. I was just wondering if there are small scripts or routines that anyone uses often in many different pipelines. Is it possible to put them on, let's say, a GitHub repo and then pull them in or to have a script stored in some common area?

(answer) Good question. So I think at the moment, most custom scripts are stored locally and executed locally. So compared to nf-core modules, which are the ship where you can download them directly. Most examples of custom scripts I think are stored locally. But with, as I mentioned very briefly just before, you can store scripts locally alongside a module with templates and that could also be downloaded at the same time as a module. But I don't think that's really been done yet. But there's nothing stopping you from copying and pasting these scripts into your own Nextflow pipeline and executing them directly from the bin directory or the templates directory. So I guess the bottom line is not that I know of, but it could be done. And there probably are examples out there that people have done this. I just don't know what they are.
(question cont.) Thanks.

[15:24](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=924)
(question) There's another question on the issue of portability when using custom scripts. Do you mean if the custom script is not within the directory tree or the base dir?

(answer) So with a custom script, when Nextflow stages all of your inputs for execution, if the script isn't defined or isn't included, it won't necessarily be found. So if you were just to include a script alongside your main Nextflow script, it wouldn't be staged because it wouldn't be able to be found when that pipeline is being executed. So if it's just in the pipeline directory, no, it won't be found. But if it's within the bin folder or the templates folder, then it can be included and will automatically be staged because it's in the path of the script or of the tool, which is a kind of... So I'm just going to read the question again. Yes, if it's not in base dir, it won't automatically be found. Although you could probably use that to specify a relative path. Although I wouldn't necessarily recommend that when you've already got the bin and templates folders available for you.

[16:39](https://www.youtube.com/watch?v=3aA5-s8PAF0&t=999)
(host) Okay. It seems there are currently no more questions. If you have more questions, you can always go to Slack and the Bytesize channel, or you can contact Chris, I guess, directly. Otherwise, I would like to thank Chris and, as usual, also the Chan Zuckerberg Initiative for funding these talks.

</details>
