---
title: 'Bytesize: nf-core modules patch'
subtitle: Phil Ewels - Seqera Labs
type: talk
start_date: '2023-03-07'
start_time: '13:00+01:00'
end_date: '2023-03-07'
end_time: '13:30+01:00'
youtube_embed: https://www.youtube.com/watch?v=7pu6Ikhi1eU
location_url:
  - https://doi.org/10.6084/m9.figshare.22231987.v1
  - https://www.youtube.com/watch?v=7pu6Ikhi1eU
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: nf-core modules patch

Have you ever wanted to tweak an nf-core module for your own use? No need to copy / paste and lose the benefits of linting and updates, instead try using the magic `nf-core modules patch` command to keep track of your modifications.

Sounds too good to be true? Join us to hear more at this week's bytesize talk with Phil Ewels ([@ewels](https://github.com/ewels))!

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:01](https://www.youtube.com/watch?v=7pu6Ikhi1eU&t=1)
Hi everyone, welcome to today's nf-core bytesize talk. My name is Phil Ewels and today I will be talking to you about the `nf-core modules patch` functionality. This is a very simple functionality, so I'm thinking today's bytesize talk will be fairly short. Many people don't know that it exists and I think it could be quite useful, especially for people using nf-core tooling and the nf-core templates for pipelines, either private or custom, which are not going to be part of the main nf-core organization. This is where this tooling really, really shines. If you want to use the nf-core templates for stuff you're doing in-house, this talk is for you. I don't have any slides or anything, it's just going to be a live demo, I'm going to walk through how I use it and try and describe what it's doing in the back end and hopefully that will make sense to you.

[1:01](https://www.youtube.com/watch?v=7pu6Ikhi1eU&t=61)
Those of you who've seen me talk before will know that I love a good live demo, things usually go wrong, but that's part of the fun of it, so let's join me on this rollercoaster. Just before I kick off, a little bit of background information, what we're talking about here. For those of you familiar with this, apologies, but just to get everyone up to the same level: with nf-core we have a pipeline template for the whole pipeline and then in the last year or two with DSL2 we've been working with modules. These are wrappers around specific tools, so this is on pipeline level and is one workflow all the way through from start to end analysis. A module is just a single tool and we have shared modules which people can collaborate on, which you can install into a pipeline. When you make changes to a module, which is a centralised module, those changes can be easily integrated into every pipeline that is using that module.

[2:03](https://www.youtube.com/watch?v=7pu6Ikhi1eU&t=123)
The classic example and the one I'm going to be working with today would be FastQC, a QC tool for DNA sequencing data used by, I don't know how many pipelines within nf-core, but very many. We've been chatting on Slack yesterday and today about some updates. There's a new version of MultiQC that's come out and it's got some new options like `--memory`, `--svg` and stuff, and we've been talking about those updates and we can just do that in one pull request, one discussion on one module and then all the pipelines can just pull in those changes across the board and get that new functionality which is fantastic. So pipelines, modules.

[2:39](https://www.youtube.com/watch?v=7pu6Ikhi1eU&t=159)
In order to make all of this system work, it's really important that the code within the pipelines, the modules within the pipeline is the same as the code in the central modules repository. That makes sense. If you want to synchronise the two, you need to keep them tightly together. What that means though is you can't change the code in your pipeline. If you do that, the nf-core code linter will start complaining and tell you that you're not allowed to do that. What people have done before is, they take the centralised nf-core module and just copy it as a local module and then they can do whatever they want with it. They can change it and the linter won't complain. The downside of doing that is you're no longer in sync, so if there are updates that come into a centralised module, you won't see them, you won't be able to pull them in and you're effectively losing that collaborative aspect which is so powerful. This is where nf-core modules patch comes in as a stopgap if you like, a way for you to make changes to central modules in your pipeline - and your pipeline alone - whilst keeping the linter happy and keeping all the functionality about updating modules and so on. Hopefully that makes sense, if you want to ask me any questions at this point, shout, otherwise we can take questions at the end.

[4:02](https://www.youtube.com/watch?v=7pu6Ikhi1eU&t=242)
Now I will dive into a screenshare and show you how this thing works. Hopefully you can now all see my setup, I'll make the zoom toolbar as small as possible. My pet pipeline is the nf-core/methylseq pipeline, it's one of the first ones I ever wrote and it's one I'm still fairly involved with the maintenance for. Hopefully everyone is familiar with the idea of the nf-core lint command which runs all the code tests on every single module in all parts of the pipeline. Today for live demo purposes I'm going to do `modules lint` which just only lints the modules and not the entire pipeline and I'm actually going to make it just the fastQC module so things work nice and quickly.

[4:54](https://www.youtube.com/watch?v=7pu6Ikhi1eU&t=294)
If I run linting, make it a bit bigger, you can see that everything's fine, my pipeline's up to date with the central dev version of methylseq and I've got a couple of warnings about this module: there's a new version of software available and there's a new version of the central module available. But they're both warnings, they're not failures, so that's my starting point. We were talking about new fastQC functionality. This is VS-Code, I'm looking at the methylseq pipeline source code now. This is not the central modules repo, this is my pipeline. If I go into modules I've got the local ones and I've got nf-core, scroll down you can see I've got the fastQC one and this is the shared fastQC module. Now I could make changes and drop into local but I'm not going to do that today. Let's say that things are moving too slowly, I want to do something here myself. What I can do is drop something custom in here, let's say I'm going to have a new input channel to handle SVGs and I want to do it just on my pipeline. I'm going to hit save. You know, assume that I'm doing some valid change here and I've tested locally, the Nextflow side of things is working fine and it's doing what I want. Now if I rerun this linting test it's going to be unhappy with me because this lint test checks the version on the web on the nf-core modules repository, looks at the code there and checks the code that I have locally and in this case it says this code does not match the version of the module that you say you have and so that's a hard failure. All continuous integration tests on GitHub will start giving a red cross and failing and this is not a good situation. This is normally where you freak out and copy it to local or something.

[6:43](https://www.youtube.com/watch?v=7pu6Ikhi1eU&t=403)
But now I'm going to do some magic. Now I do `nf-core modules patch`. I run this command, it asks me which module I want to do it for, the FastQC, that's where I've made my changes and it just very quickly spits out some content. First things first, this is a diff, so this is where it's looked at the remote file and the local file for any changes. You can see it's picked up here that there's some code in my local chain copy which has changed. This looks right, this is what I just added, so it says there's an extra line here. Now these diff files are really cool because with diff you can generate these diff files or patch files and you can apply them on top of other files, so we can reapply this change at any time. That's what we do, we save this diff and if I go to git status you can see I've got changes to my modules.json file which is used by nf-core to track the synchronisation between my local pipeline and the shared modules repository. I've got the changes in the FastQC file which is the thing I just edited and saved and I've got a new file here called Fastqc.diff. If I go back into the VS-Code we can see that this diff file is just what was printed to the console here and it's just saying alongside the FastQC module, hey I've got some local changes here. Then if I go into modules.schema you can see if I find FastQC that we've got a new line that's been added in the JSON file here and it's just telling nf-core that there is a patch file that exists in this location.

[8:19](https://www.youtube.com/watch?v=7pu6Ikhi1eU&t=499)
Okay so great, what does that do? I can add all of this, let's make a new branch. Now if I go back to the lint commands, `nf-core modules lint`, which was failing, we're back to our starting position. Everything's fine, everything's happy. Now just to explain what's going on here, in the background I've still got those local changes but in the background when I do `nf-core lint`, the nf-core code is fetching the remote version from the nf-core modules repository but then it knows I have a patch file, that diff file. It stores that local copy it's got from the web, it applies the patch file on top and then it compares. That's why there are no changes. If I make some more changes in here again, so `val foo`, then that's not going to be in the patch file and it's going to fail. In fact it did a hard failure where it couldn't even figure out what was going on if I do it in a different place here. Then it will just fail again and say that something has changed. Then I could run `nf-core modules patch` again, it will update that diff file but now there's new changes are covered by the diff and everything will work again. Hopefully that makes sense.

[9:50](https://www.youtube.com/watch?v=7pu6Ikhi1eU&t=590)
What's cool is it's not just linting, this applies to. The same process works when I update modules. We've got this lint warning that there's a new version of this module available. I can also do `nf-core modules update`, let's just do FastQC and hopefully, yes, there we go. It has updated the module for me, so it's gone to the nf-core remote with a shared one and it's updated my local copy and then it's managed to still reapply the patch file on top of a new updated version and save that there. `nf-core update` still works, which is like magic, I think. We can see these are all the new changes that have come in. There's a new git_sha for this module for the latest version when you can see the changes which happened when I updated the module.

[10:49](https://www.youtube.com/watch?v=7pu6Ikhi1eU&t=649)
Right, hopefully that's all you need to know, everything works beautifully, but I thought for completeness, I would also show you one small complication of when things could go wrong. We got a hint of it a second ago, actually. Something that could happen is if I do reset, so just go back to when we first made the patch before we updated. Now I can add a different change here. Now I'm going to add `--svg` onto the FastQC command itself. I execute `git diff`and you can see it's the same, `nf-core modules lint fastqc`, you can see it's the same, this is all the same, `nf-core modules patch fastqc`, yes, regenerate the patch, okay, so now our patch file has got two changes here, that's good, and `nf-core modules lint fastqc`. It's the second change, but everything is exactly the same so far. Now the tricky bit comes now, if I go to `nf-core modules update fastqc`, just like I did before, it will fail.

[12:14](https://www.youtube.com/watch?v=7pu6Ikhi1eU&t=734)
Now what's happened here, I've got a couple of warning messages saying it's failed to apply the patch. You will have to apply the patch manually. This is a little bit like when you are working on code and you put in a pull request and you get a merge conflict, there have been changes that have happened on the central nf-core modules repository, and there have been my local changes which I've done with those patch files, and the tool can't automatically figure out how to reconcile those two changes. What it's done is it's just clobbered my local changes. If I go in here, you can see it's made the updates, but I've lost all my custom changes and it's just overwritten it with what was on the remote copy. All I have to do is I have to go back in and I have to just recreate this patch. I can go back in and go `--svg` and `val svg`, and then rerun nf-core modules patch. That's fine, so just bear that in mind. Sometimes when you do updates and you have local patches, you might need to do a little bit of fiddling. Just be careful with always running git commit before you do stuff, because then you don't lose anything and you can easily see which changes are happening.

[13:28](https://www.youtube.com/watch?v=7pu6Ikhi1eU&t=808)
Right, that's my live demo. Hopefully everybody followed along with that and it made sense. Yes, sorry. Fran says in the comments that when I said earlier about shouting, you can shout literally if you want to, but also then just you can ask. Happy to take any questions and hopefully this will be useful for some of you.

(host) You can unmute yourself now if you want to ask questions. Lots of happy people, very few questions.

(speaker) That's good. It's quite a nice, neat little small topic to discuss, so I had plenty of time for it.

[14:20](https://www.youtube.com/watch?v=7pu6Ikhi1eU&t=860)
(question) Sure, but what happens if you realize that you actually have more changes than you actually want to include, like if there's something that make the whole pipeline fail or whatever. Can you undo the patch.

(speaker) When you do the update, you mean?

(question cont.) No, when you have written something, you get errors and it's like, ah, but this is because I did something manual that I want to be different. You do a patch and then you realize, oops, actually that was something completely different and messed everything up.

(answer) Then it's no difference if you made changes in any other way. You just go back and you look at your git changes and you revert in git or whatever. This part of the workflow is specific to just the nf-core tooling, just the linting, the updating. It's coming in at the end once you've already fiddled around and made Nextflow in your pipeline work properly in the way you want it to.

(host) Cool. Then thank you very much. Also thank you everyone for listening and as usual, I would like to thank the Chan Zuckerberg Initiative for funding our bytesize talks and "Hello" also to Maxime.

</details>
