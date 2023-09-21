---
title: 'Bytesize 12: Template sync - how to merge automated PRs'
subtitle: Phil Ewels - SciLifeLab, Sweden
type: talk
start_date: '2021-05-11'
start_time: '13:00+02:00'
end_date: '2021-05-11'
end_time: '13:30+02:00'
youtube_embed: https://youtu.be/-CZKoo5Y_08
location_url:
  - https://youtu.be/-CZKoo5Y_08
  - https://www.bilibili.com/video/BV1D64y127yJ/
  - https://doi.org/10.6084/m9.figshare.14572866.v1
---

# nf-core/bytesize

Join us for an episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation.
Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 12: Template sync - how to merge automated PRs

This week, Phil Ewels ([@ewels](http://github.com/ewels/)) will present: _**Template sync - how to merge automated PRs**_

This will cover:

- Introduction to the nf-core pipeline template
- Overview of the `TEMPLATE` branch and sync theory
- What the automated `nf-core sync` does
- How to merge sync PRs

The talk will be presented on Zoom and live-streamed on YouTube:

- YouTube: <https://youtu.be/-CZKoo5Y_08>
- Bilibili: <https://www.bilibili.com/video/BV1D64y127yJ/>
- FigShare: <https://doi.org/10.6084/m9.figshare.14572866.v1>

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:46](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=46) Hi everyone, I’ll be going over how template synchronization works. This talk is targeted at pipeline developers within the nf-core ecosystem and for those who are interested in building their own pipeline based on the nf-core toolsets.

[1:03](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=63) First, a refresher here. All of this is based very heavily on Git and GitHub, so let’s very quickly go over some of what we’ve covered before [Bytesize#4](https://nf-co.re/events/2021/bytesize-4-github-contribution-basics). Git is a version control code system that has a bit of jargon that you need to know. You do `commit`, which are little packages of work, so you do some edits and then bundle them together and then say `commit` and that’s like checkpoints in the timeline of your code. They go in a kind of straight line, and then you can make branches, where you can work on different things at the same time.

[1:44](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=104) You can carry on working in that branch.

[1:47](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=108) And when you want to, you can merge that branch back into the main line of the commit history.

[1:55](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=115) You can fork an entire repository, which basically copies a repository of code into your account. So for example, I have forked an nf-core pipeline to my personal Github account.

[2:08](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=128) And you can carry on doing development on your fork. When you want to contribute what you’ve done back to the main repository within nf-core, you open up a pull request that is then merged into the main line of commits. Then everything you’ve worked on becomes a part of the main repository.

[2:29](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=150) All our nf-core pipelines have at least three branches called `dev`, `master` and `TEMPLATE`. The `master` branch is for default (after the first release), and should only ever have the latest stable release code. That is so that if someone runs Nextflow with the pipeline name, but does not specify a release, it can grab the default branch, which has the latest stable release.

[3:00](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=180) The `dev` branch is where we do all the development work, so things change there all the time. When we’re ready for a release, we merge `dev` into `master`.

[3:14](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=194) Finally, we have this kind of special branch called `TEMPLATE`, and that’s what I’m going to be talking about.

[3:24](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=204) So this is one of the figures from the [nf-core paper](https://rdcu.be/clDKw), and it describes how automated sync works from nf-core `TEMPLATE`. One of the main things we offer within the nf-core organisation is this pipeline template. You can do `nf-core create`, where it asks you a couple of questions i.e. a name, a description, an author, and it creates a minimal but functional pipeline. It has a lot of the boilerplate code that we think is best practice for a Nextflow pipeline. So it sets up lots of config stuff and it creates all the code structure into different files, it sets up automated tests and Github actions, continuous integration etc. But that `TEMPLATE` never stands still, we’re adding to it all the time, and it grows with the community and experience. So if we fix something in one of the pipelines, which is general, we try and bring that back into the main `TEMPLATE`. But of course once you’ve done `nf-core create`, your repository is then separate and no longer tracks along with the nf-core template, and that’s what we’re trying to bring in here. So the trick we’ve figured out is that when we make your pipeline, you straightaway make your commits where you haven’t done any changes - it’s just a template. You have that initial `commit` with the `TEMPLATE`, and we branch at that point to a new branch called `TEMPLATE`, which has the same shared `commit`. You can then just forget about that `TEMPLATE` branch and carry on doing whatever you’re doing in `master` or `dev`, but when there’s a new release of the nf-core `TEMPLATE`, we create a new `commit` along this special `TEMPLATE` branch. Now we’ve got one `commit` in this to commit, and the only changes in this branch are the things we updated in the main `TEMPLATE`. This is always the case for the `TEMPLATE` branch, it only ever has the vanilla template from the main nf-core `TEMPLATE`, and the commits are only ever just updates from there as we do new releases of tools. When you get a new `commit`, then we can do a pull request from that, `branch` it into the `dev` branch so that anything that changed in the boilerplate code is brought into `dev`. But it never goes the other way, dev is never brought into `TEMPLATE`.

[5:39](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=337) So hopefully that made sense.

[5:44](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=344) Finally, we automate this process, so that when we do a new release of the tools repository, all these changes can happen automatically and the pull requests are automatically opened.

[5:57](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=357) So what does this synchronisation process do step by step? It’s pretty simple. First of all, it gets a copy of your pipeline and checks out the special `TEMPLATE` branch and deletes everything in that pipeline. Then it creates a new pipeline from scratch using `nf core create`. This time of course what you’re using is the latest release of the nf-core tools package, and so it uses the latest version of a `TEMPLATE`. So this new version now has all the modifications that we’ve made since the last release. Those changes, if there are any, are committed into a `TEMPLATE` branch (we’ve recently added a feature where we make a special branch and I’ll go into why). Then we open a `pull request` with these changes to the dev branch, and Git takes care of the rest. You’ve been doing thousands of `commit`s on the `dev` and `TEMPLATE` branches, and a whole lot of other things, but Git picks what’s relevant. So then you review that pull request and have the latest changes.

[7:33](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=453) We’ve got something special for you today, which is the `nf-core/tools` here. We’ve been working hard on a new release, and in fact it’s all prepped and ready to go.

[7:46](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=466) Here we are on the release page. I’ve just written the changelog and we’ve merged the master, all tests are passed, so I’m going to hit release now. We will see if the synchronisation process works (automated sync doesn’t work completely flawlessly all the way and will mostly fail).

[8:09](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=489) So I hit publish, and there you go, tools is released!

[8:13](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=494) Now if I go to the Github actions page, you can see lots of stuff is automatically happening here because of the release. The one that we are interested in is called `sync template`.

[8:24](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=504) Here we go, the sync has started. The first step is to find a list of all the nf-core pipelines which it gets from the [nf-core website](https://nf-co.re/pipelines), then kick off jobs.

[8:34](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=514) There’s a synchronisation job for every single one of the nf-core pipelines. It goes through all of these and synchronises them and opens up a pull request for every single one.

[8:47](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=527) There are now so many, that Github thinks that we are abusing the API, so it stops all our pull request, and that’s what I’ve added code to fix.

[9:05](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=545) We will see how many red crosses we get.. You see that it’s going through all the steps I just talked about. It has checked out the code, installed nf-core, then Nextflow, running the synchronisation, got all the variables it needs from a dev branch, checked out the `TEMPLATE` branch… and then crashed horribly.

[9:32](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=572) But you can see it committed changes and pushed to the `TEMPLATE` branch on the remote (found a new bug). It made a new branch, pushed it to remote and then tried to make an API and then found something wrong in here (a bug), and it’s crashed now, but it worked before it crashed.

[10:10](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=610) So thankfully, we have our automated pull request. Now if I go to this repository, we can have a look at the branches. You see the template branch here, and see a `commit` for each one of the releases of tools.

[10:25](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=625) This is the new commit that we just pushed, and these are all the things that have changed in the template since the last time.

[10:33](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=633) Now if I go to a special branch called `nf-core-template-merge`, it looks exactly the same and that’s where the pull request comes from.

[10:42](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=642) So this is where you start if you’re a pipeline developer. You will just have got a pull request in your inbox that you will go in and see something like this that hopefully explains what to do.

[10:54](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=656) Now, there’s some stuff running here and this is where it gets a bit more interesting nearly all the time.

[11:02](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=662) There are some merge conflicts because you’ve been working on all the pipeline files. Then in this parallel branch, there are changes to the same pipeline files. So Git does well in figuring out what should go where, but it does usually need a little bit of help, sometimes it’s just one or two files or it can be half the pipeline. So there’s almost always some merge conflicts to fix.

[11:32](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=692) There are a couple of ways to do this if they are small and simple.

[11:37](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=697) This button here on Github will not be greyed out, you can click on that, which will open a browser. I can look if I can find where it is allowing me to do this.

[11:56](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=716) There you go, we’ve got some conflicts here, but they are a bit simpler and this button is not greyed out. I can click on that..

[12:01](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=721) ..and it will take me to the files which have got the merge conflicts. You can see it highlights the merge conflicts. So it says that this line has probably been edited in both.

[12:15](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=735) This looks more up-to-date and more like what we want to keep, so I delete this part of the merge conflict, and delete that marker and mark it as resolved. Then I just go on to the next one, and work through the rest. This is only really viable if the merge conflicts are few and simple because you need to manually deal with the merge conflict markers, and it’s quite easy to make a mistake. So generally, I would avoid doing this very often.

[12:44](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=764) Very often that button will be greyed out if there’s too much for Git to do on the web interface, and you will need to use the command line instead. You can basically do that on a branch, but I would recommend doing this on your fork instead.

[12:58](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=776) So the process would be to first update your personal fork at the repository against the `dev`, make a new branch to do this work, and then pull the changes from the upstream repository on the `TEMPLATE` branch.

[13:17](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=797) If you watched [bytesize#11](https://nf-co.re/events/2021/bytesize-11-dev-envs-workflows), I talked a little about how to pull in changes from upstream repositories. So check that out if you’re not quite sure. But typically, you have two forks, you have origin which is your fork, and upstream which is the nf-core fork. So you can do ´git pull upstream template´ and that will pull the changes from the `TEMPLATE` branch (there will be a bunch of merge conflicts but you can fix those). You can commit them, push to your fork, and then open a pull request from your fork into the main repository.

[14:03](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=841) Remember that `pull request`s are updated dynamically, so when you do that, you will have pulled this `commit` into your personal fork, and when you open your pull request, it will have the same `commit`, and your `merge conflict` fixes. Merging that pull request will result in it automatically showing up as merged because git will see that this `commit` was to the `dev` branch.

[14:27](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=867) So basically if you’re going to do that, don’t worry about this pull request. Don’t close it, just leave it there because if you run out of time or don’t finish your personal fork, it’s clear that there’s a template update waiting to be done.

[14:50](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=890) Most importantly, follow the documentation.

[14:55](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=895) If I go to [developers template synchronisation](https://nf-co.re/docs/contributing/sync), there is a step-by-step walkthrough on how to do everything that I’ve just described.

[15:09](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=909) It says how it works, and then we say how to actually fix the minor conflicts and how to deal with the things I just described during my talk. There are example commands etc, so I follow this documentation myself as well.

[15:29](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=929) It also has some other stuff about how to do much more difficult types of synchronisation, and a few more things on how to fix stuff when they go wrong. Hopefully most of that is behind us now that we have this new system. The reason this comes from a branch is that if you fix the merge conflicts using this little button. This will pull all the dev history into this branch as part of fixing the major conflicts. If these come from a `TEMPLATE` branch, that’s bad because `TEMPLATE` is where your development code is, and the next time you run a `template sync`, it will end up deleting the whole pipeline and creating a template. Git then will think that it is what you wanted to do, but it isn’t. So that is why we want that `TEMPLATE` branch to be completely clean, and have started using these branches although that are discarded afterwards.

[16:33](https://youtu.be/-CZKoo5Y_08?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=993) One of the things that came out in the last version of nf-core lint for pipeline linting is an option called `--fix` where some pipeline tests. For example the important ones which check that a file hasn’t been modified from a template, so if there’s some core file which you never need to touch, it will compare the two, and if you do `--fix` then the test name will delete that file and create what it expects. Make sure to check the `git diff` after you run this, mainly to check that it hasn’t deleted any of your custom code. This can be much quicker than trying to resolve merge conflicts, you can just run `nf core lint` and you will end up with the right file. It’s just always nice to go through and delete any of the branches that you don’t need anymore.

</details>
