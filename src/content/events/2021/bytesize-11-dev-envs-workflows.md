---
title: 'Bytesize 11: Development environments & workflows'
subtitle: Phil Ewels - SciLifeLab, Sweden
type: talk
start_date: '2021-05-04'
start_time: '13:00 CEST'
end_date: '2021-05-04'
end_time: '13:30 CEST'
youtube_embed: https://youtu.be/XB96efweCLI
location_url:
  - https://youtu.be/XB96efweCLI
---

# nf-core/bytesize

Join us for an episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation.
Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 11: Development environments & workflows

This week, Phil Ewels ([@ewels](http://github.com/ewels/)) will present: _**Development environments & workflows**_

This will cover:

- Git + GitHub configuration
  - Workflow for forking main repos
  - Configuring git (remotes / pulling upstream updates)
  - Convenience command line aliases + functions
- VSCode configuration
  - Code linting (Black / Prettier / EditorConfig / Markdownlint)
  - Handling merge conflicts in VSCode
  - Other useful VSCode plugins
- Pipeline testing
  - Running tests locally / reading the CI test logs
  - Ignoring `nf-core lint` tests
  - Troubleshooting software builds (conda / building Docker images locally)
  - Troubleshooting inside containers

The talk will be presented on Zoom and live-streamed on YouTube:

- YouTube: <https://youtu.be/XB96efweCLI>

You can see the slides on HackMD: <https://hackmd.io/@nf-core/bytesize-11#/> (shown below)

<div class="ratio ratio-16x9 border shadow">
  <iframe  src="https://hackmd.io/@nf-core/bytesize-11#/" allowfullscreen></iframe>
</div>

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[1:04](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=64) So today’s talk is in response to the suggestions for talks for this series that we received from you. Some of you said that you would find it interesting to see how one or more of us who set up things on nf-core organise ourselves, with all the tips and tricks we have to make our day-to-day working easier. I have tried to build this around a talk that Alex Peltzer gave a few sessions ago (see [Bytesize#4](https://nf-co.re/events/2021/bytesize-4-github-contribution-basics)). If you’re completely new to working with GitHub, please check out his talk first because he covers a lot of jargon that working with GitHub entails. I’m going to cover slightly more advanced things during this talk. This is based on my personal set up, and there are ways to discover and figure things out as you go along. But hopefully, some of what I show here will be helpful. I realise that the term workflows is probably a bit misleading in the context of nf-core. Workflows here is meant to be a demonstration of how I work, how I get into a project, start working with a Git repository etc.

[2:52](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=172) I will cover three sections and I’m going to start by talking a little about GitHub configuration: how I start on a new project and how I work collaboratively on repositories. The one thing that comes up a lot on Slack is how to pull in updates other people are doing on the `dev` branch, and this is relevant when you have lots of different people working on the same code base, especially if you’re doing this for the first time. Then I’ll briefly touch upon my code editor; I used to use [Atom](https://atom.io/) for a long time, and have recently switched to [VSCode](https://code.visualstudio.com/). Finally, I’ll cover a little about pipeline testing and how I do some debugging.

[4:02](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=242) There’s going to be a mix of slides and live demo here. So Git and Github configuration: four topics here. Firstly, getting started with working on nf-core or any repository that you haven’t created yourself to start with, but you’d like to contribute to. This is the case for most nf-core pipelines, most of the time. So, you’re using a pipeline, and you’ve found a bug or want to add a new feature. Let’s talk about how to get that pipeline off the web and onto your system, synchronise the updates, and contribute to other people’s pull requests with a live demo.

[4:46](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=286) Just before we started, I dug out an nf-core pipeline that I have never worked on before. `nf-core/cageseq` is one of those pipelines that’s one of the newer pipelines and I’ve never cloned it, or contributed to it.

[5:09](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=309) The first thing I want to do is to make my own fork from the head nf-core repository.

[5:24](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=324) I do that by clicking `fork`. It asks where I would like to have it forked to, and I’d like to have it forked to my personal account. So I do a couple of things with this new repository, and the first thing I do is to drop the nf-core prefix. So I always do this mainly for myself so that it’s clear to me what it is. The other thing I do is to grab the URL of the main one, and update this “about” field and say main pipeline here, and put the URL for the main repository here. Then I uncheck all this stuff, turning off as many features as possible because I don’t want it to collect issues, wikis, projects, or anything like that. So I turn all of that off.

[6:16](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=376) I’d like to make it clear as possible that this is not the main repository and doing this stuff doesn’t affect the upstream fork.

[6:28](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=388) So the next thing I want to do is to work on this code locally. I do that by cloning the repository. So I do a little drop down and copy this onto the clipboard. You can use either `https` or `ssh` (I personally always use ssh) - you need to set up GitHub for that, but it’s not too difficult to do, and I will show you an easy way to do that.

[6:51](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=411) So I copy this URL and then I go to my GitHub, do `git clone` and then I paste the URL. Then we get that into a directory. I usually rename it.

[7:15](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=435) Then if we just go into this directory, I’ve got the files there. I can do `git log` to look at all the recent changes, `git status` to check what I’ve done, and `git branch` to see all the branches and so on.

[7:32](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=452) So I can push my changes to my fork of the repository. Now let’s assume I want to work on something here, so I’ll make a new file.

[7:47](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=467) Now if I do `git status`, I see that there’s a new file; `git commit -a` and because it’s a new file, I need to do `git add text.txt`. Then I do `git status`(you see it’s staged here), but I need to commit it, so I do `git commit -m`, m for message “This is a new file”. I can do `git push`, and it will push that file to my fork of the repository.

[8:15](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=496) The main nf-core fork still doesn’t know anything about it. So to start doing things in the main fork now, I can open up a new `pull request`.

[8:25](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=505) I’ve made an easy mistake here, which is that I was working off the `master` branch - the default branch, and so the code I was changing was quite well behind the `dev` branch, where the latest version of the code is going to be. Sometimes that doesn’t matter, but sometimes it could be a bit of an issue. In this particular case, I’ve just created a new code from scratch, but if I was editing code that had already been changed on `dev`, it could be a mess.

[8:48](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=528) So what I needed to do was to pull the changes from `dev`. If I just undo this by `git reset --hard HEAD -1`, that just takes my repository back one commit and deletes everything. Then I do `git push --force`, which overwrites the remote.

[9:14](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=553) Now I want to pull in the changes from the main nf-core workflow, so I copy this URL again. What we want to do is we want to add a second remote. So when you clone a repository, it sets up a remote for the web repository at github.com, and you can see that if you do`git remote -v` (v for verbose), you see that I have a remote called origin and the URL that I pasted in. I want to tell my local copy of the code where the upstream version is, I do that through `git remote add`. I call it upstream here, but you can call it whatever you want (convention is usually upstream). I’m going to paste the nf-core URL. Now if I do `git remote -v`, you can see that I’ve got both sets.

[10:01](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=601) Now I can do `git pull` and this gets the name of the branch I’m on. I can do `git pull upstream` and the name of the branch (I’m interested in getting stuff from `dev`) to pull up the new changes from the `dev` branch.

[10:25](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=625) `git push` pushed that to my fork, which is now up-to-date with the `dev` branch instead of with the `master` branch. I can go ahead and do changes and then make a `pull request`. It is important to pull changes like this every time you come back to a repository because it’s a good way to avoid merge conflicts and so it’s really important to get used to doing this. It is a key concept when working collaboratively and when a lot of changes are being pushed all the time.

[11:22](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=682) You can also do this if you’re working on a branch, you can still pull from whatever upstream branch you want.

[11:28](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=688) One tip - something I’ve reinvented myself - because I do those sequencer commands so frequently. I’ve written myself a little bash function, which I call g update, so I do `gupdate dev`, each time I start working on a project on the terminal. It also prunes branches and does a couple of other things. The [slides](https://hackmd.io/@nf-core/bytesize-11#/) for this presentation have links to where I keep my snipped codes, so do check them out if you’d like to use them. You can also make your own version of the same thing.

[12:04](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=724) One final thing with Git workflows is again to avoid merge conflicts. It’s nice to work on branches because it allows you to do a packet of work on a pipeline, make a `pull request` that can go back and forth during code review and that could take a while, but you might also want to do another fraction of work on the same repository at the same time. If it’s all together on the same branch, it can be very difficult to execute. So instead, it’s a good idea to start out on a new branch using `git checkout -b`, which is a shortcut for creating a branch and checking out to it.

[13:03](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=783) I almost never do any work on `master`. Now I can do my normal work, so `echo “Test” > test.txt`, `git add test.txt`, `git commit -m My new file`. When I push this, it will complain because on my fork on GitHub, I don’t have a branch called “My new feature”, and I haven’t set up where to push that. So the first time you try and push a new feature to a new branch, you might not be able to. In that case, you just copy that command and paste it. Alternatively, you do `git push --u origin` (origin is the name of my remote for my fork). So that has now created a branch on my fork with my new feature.

[13:56](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=836) So there it is.. You can see that it’s got the test.txt file and what’s especially nice is that GitHub actually even responds with a link to make a pull request. I can just hold the command key (`command`) on a Mac or the control key (`ctrl`) on a PC, and click to open up a new pull request straight out of the terminal. So this is the main workflow that I follow every day for working on code-base, making pull requests, updates etc.

[14:30](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=870) It’s taken me a long time to get to a workflow that I really like when working on other people’s pull requests. So for example, if I go to nf-core/tools and look at these pull requests, I notice that there is a pull request here that someone has made. It has a bunch of merge conflicts and some small bugs that would take longer to write comments on than to fix yourself. It could also be that you’re working with someone else on a pull request but they’ve initiated it from their fork. There are lots of different ways or reasons why you might want to push your code to someone else’s fork. This relies then on the person who created the pull request to have a box (“allow maintainers to push to this pull request”) checked, but by default it’s checked for personal accounts.

[15:34](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=934) What I use is a tool released by GitHub called `github cli`, which you can install on any system in a variety of ways. This gives me a new command called `gh`, and it’s both very powerful, and flexible. So in this case, I’m going to go to my tools repository, then `git checkout master` to go back into the `master` branch, and then I’m going to do `gupdate dev`, so that it pulls out a bunch of stuff that’s been merged into the `main` branch. Then I’m going to copy the pull request number here, and then I’ll do `gh pr checkout` because I want to check that code out on my local system. I’ve already set up the CLI with authentication and it knows that the nf-core repository is like the main one where I’m interested in tracking pull requests. So it just knows what to do, it just gets everything for me and pulls this code out into a new branch that I now have locally and can work on. I can even do `git commit` `git push`, and it will turn up on the other person’s repository on their fork and feature in this pull request. Now multiple people can work on the same progress space as long as they have the right access. This is incredibly powerful.

[17:10](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1030) I’ll give you a quick example. So here we have a merge conflict, a classic changelog where two lines have been added by two different people in parallel and it doesn’t know which one to use. But we can see we want to keep both, so I could just do this on GitHub on the web browser by just deleting these lines and then clicking `mark as resolved`. Now for the sake of providing a demo here, I’m going to do it locally. This is also because for some merge conflicts, it can be more complex to merge on the web interface.

[17:58](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1078) Use a code editor! So what I’m going to do is to pull in the changes from my `master` branch, which should generate the same merge conflict. Then I’m going to try and fix that. I’ve got VSCode installed, so I’m going to do `code .` to open this directory, and then if I click the git tab, I can see merge changes (with an exclamation mark). This is the file that needed my attention. I can click on that icon, it takes me to that file, I see we have a merge conflict and it’s highlighted, so I can scroll through this file and see all the merge conflicts. What’s nice is that there are buttons along the top saying accept current, which is green, accept incoming, which is the one I’m trying to merge in from `master` and is highlighted in blue, and then there’s accept both, which is for more complicated things. I want to keep both, so I can just click and it’s done.

[19:10](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1150) Then I go do `git status` - I use another shortcut here - `gs`, which is more succinct. You see it’s all fine, then I do `git commit`, `fix merge conflict`, `git push` even though this directory is my fork. This pull request is basically the other person’s pull request, but my commit ends up on this pull request.

[19:51](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1191) Now when I’m done I just go `git checkout master`, and if I want to I can `git branch -D` because it’s not yet merged. I can get ride of that branch which was pulled up.

[20:09](https://youtu.be/XB96efweCLI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1209) So that was the main part of the talk. I’ve covered Git, GitHub, how I work with other people’s code, remote code, VSCode configuration, merge conflicts etc.

</details>
