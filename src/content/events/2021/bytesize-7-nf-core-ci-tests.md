---
title: 'Bytesize 7: Making the CI tests pass'
subtitle: Phil Ewels - SciLifeLab, Sweden
type: talk
startDate: '2021-03-16'
startTime: '13:00+01:00'
endDate: '2021-03-16'
endTime: '13:30+01:00'
youtubeEmbed: https://youtu.be/U9LG_mMQFMY
locations:
  - name: Online
    links:
      - https://youtu.be/U9LG_mMQFMY
      - https://www.bilibili.com/video/BV155411P7ES
---

This week, Phil Ewels ([@ewels](http://github.com/ewels/)) will present: _**Making the CI tests pass.**_
This will cover:

- How automated tests work
- nf-core lint
- Code syntax / formatting
- Editor integration

The talk will be presented on Zoom and live-streamed on YouTube:

- Zoom: <https://zoom.us/j/95310380847>
- YouTube: <https://youtu.be/U9LG_mMQFMY>

<details markdown="1"><summary>Video transcription</summary>

:::note
This text has been edited to make it suitable for reading
:::

The talk today is going to be less slides and more focused on a live demonstration because it’s probably best to learn this by doing. So please bear with me.

So a quick introduction. What is CI? CI stands for continuous integration, and basically that’s just another way of saying automated tests.

So every time we push a change to the code for any repository that has it set up, push a change to Github, do a `git commit` or a `git push`, that triggers some sort of CI. This will typically result in a test being run and a report that indicates whether that test has passed or failed.

[1:46](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=109) The test itself is run by a dedicated service; the ones that have been most popular in the open source community are ones called Travis CI, which we used to use years ago on nf-core, Circle CI, which is used quite a lot and I think Sarek and BioConda use it, and then there’s the one we use at nf-core called GitHub Actions.

[2:46](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=166) GitHub Actions is an added feature of GitHub itself, as an opt-in. Now these are super-flexible and powerful. It’s not just tests that you can run, but all kinds of different things. It doesn’t have to be code changes, you can have events such as opening a pull request, making a release, or even clicking a button on Github webpage if that’s configured.

So lots of different inputs to trigger things and we can have lots of different things that result. This is what I’ll cover during this talk today.

[3:12](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=193) So you’ve likely often heard about continuous integration as CI/CD. CD is continuous deployment, which is the second part.

[3:34](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=215) I’m just going to go into the live demo directory. I’m hoping it’s going to go better than the Tesla demo in the gif on the right here.

[3:52](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=232) So, here is a pipeline. This is an nf-core pipeline. We have a special one actually, called the test pipeline. For all purposes, it behaves exactly like every other nf-core pipeline, so we can use it as a test-bed to check the code.

[4:11](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=251) I should note that much of what I’m going to demonstrate here revolves around the nf-core tools helper package and we have a big release stacked up and ready to go.

[4:34](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=274) As a a result, the demo I will be showing you will be with version 1.13 dev, which is about to become version 1.13 stable. That’s why things are a bit funny here [4:43](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=283).

[4:46](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=286) Then if you try and replicate what I’m doing right now with the production, with the main stable tag, it won’t work yet. But hopefully when we release this.

[4:52](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=292) So anyway I’ve got my nf-core/pipeline here, and I’ve made some changes such as updating the `CHANGELOG` here [5:00](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=301).

[5:06](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=306) I’ve forked this pipeline, so you can see this is my Github username, and this is this pipeline that you can see is forked from the main one.

[5:13](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=312) I’ve actually just pushed my change before the talk started, so if I go to the list of commits [5:16](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=318), you will see the latest commit that I pushed here [5:21](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=321).

[5:23](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=323) If I remove a bunch of to-do things, this is pushed to GitHub [5:30](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=330), and you will see that next to each one of the commits here in this log is a little tick or a cross. These are a summary of all the different tests. I can click on this little tick [5:42](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=342) wherever I see it, and I’ll see a list of the different tests.

[5:46](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=346) You see here that I’m being told that some of the CI tests are not successful.

Now if I jump here [5:51](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=351), I can actually go through to this panel [5:53](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=353), which is a part of Github actions and it tells me all the different tests which are running.

Now I can see a breakdown of all the different tests [5:58](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=358), and whether they worked or not.

[6:03](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=363) The one that’s jumped to this Markdown one, and this is a part of a group of nf-core linting.

What this does is that it looks over markdown code, which is written and we validate it against a code formatting check.

This is something we added a little while ago, and it doesn’t really matter how you do your markdown. In practical terms, it renders the same when people read it.

But it’s really helpful for us to code in a standardised way because we have so many contributors with nf-core.

It allows us to follow some sort of common guideline, and so we make heavy use of what we call code linters or code formatters, and this is what we are running into here [6:42](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=402).

If you’re new to nf-core, it’s likely that you’ll hit this if you haven’t used it before. You’ll just add something to a `CHANGELOG` and it’ll look fine to you. Then you’ll push it and see a test failing.

So what this test does is it runs a command called ´markdownlint´ [6:56](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=416), which is a package on [npm](https://www.npmjs.com/package/markdownlint) which you can install yourself.

[7:00](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=420) You can see all the different stages of this and the bit that’s failed its run. I can also see that the command has run.

[7:06](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=426) So I can run this locally if I want to, and it should give me exactly the same output.

[7:11](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=431) And you can see that it says on this file `CHANGELOG` on line 6, that there’s something wrong with this rule in this linting tool where it says headings should be surrounded by blank lines.

[7:28](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=448) You can write markdown or .yml or these other linters that we have and just see if it fails and then fix it manually, but that’s kind of dull.

[7:43](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=463) So what we recommend is to set these up on the command line yourself, and all these tools have options to fix these things in place. But even better is to install plugins, which are a part of your browser.

[7:52](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=473) So I have a plugin for vscode, which is what I am using here (and is fair for markdown linting).

[7:58](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=478) In fact that is why I have this wiggly line along here because it says there’s something wrong with it. But best of all, it has automatic formatting built-in.

So if I save this file [8:06](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=486), it fixes it for me automatically. I didn’t need to do anything.

So once you have your environment set up like this [8:10](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=490), you just kind of forget about it and it always works because every time you hit save, it automatically fixes everything for you.

So I recommend setting this up. We’ve recently done set this extension pack up (with included extensions), and one of them is called `markdownlint`.

[8:35](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=515) So if you just type in nf-core into the vscode package, you find this [8:39](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=519).

Install it and activate all the things. One of them will be markdown, and you’ll get this kind of magic behaviour.

[8:46](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=526) OK, so that’s fixed. If I go back and say `git status` [8:52](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=532), and I can do `git diff CHANGELOG.md` [8:56](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=536), and there you go, there’s my extra white line [8:58](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=538),`fix markdownlint`, let’s push that change and we’ll see if that makes everything happy.

[9:09](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=549) Go back to my pipeline, and now I can click on actions here [9:12](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=552), and I’ll just see everything that runs.

You can see that this is a new action named after the commit running here [9:16](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=556), or I can go to the little status icon next to the commit itself, and you can see this is still running at the moment, but I can click on it and see it running in real time.

[9:30](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=570) So it’s just installed markdown lint and it has run markdown lint, and there is no output, which is a good thing.

[9:35](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=575) And now I have a green tick next to markdown, so we have solved problem number one.

[9:43](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=583) So the first class of continuous integration tests that we have is code linting or markdown yml, they’re all the same. Set up in your code editor and forget about it.

[9:54](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=594) If you’re writing with nf-core/tools, there’s one for python as well, which is called [black](https://pypi.org/project/black/), but they’ll do the same thing.

[9:54](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=598) Now there’s some other tests to be run, which are not quite as trivial, and the most important one is the nf-core one, which is the one we’ve built ourselves for the community.

[10:10](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=610) This has been around since the start of nf-core and what we realised was that it was almost impossible to manually check that everyone was adhering to all the guidelines and best practices exactly as we wanted.

It was just too easy to miss stuff. So we like automation and have built a tool that checks the code in your pipeline. We update this continuously and that means that every time there’s a new release, it will come up with new tests. So pipeline developers, your CI test might fail until you’ve updated your pipeline.

[10:44](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=644) That’s what we want. So here, you can see an example of a pipeline that has run this test and failed. So I can click on it and we can see the results now.

[10:55](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=655) This is a command line tool and you can also run this locally.

[11:02](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=662) So if I go to the pipeline and do `nf-core lint` and give it the pipeline directory, in this case a dot, it runs all the checks and does a few things about Conda packages and things and figures that out.

[11:21](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=681) We get the same results here.

[11:23](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=683) I recommend running this locally, just to check. One of the really nice things is that the output here has interactive hyperlinks built into it.

[11:34](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=694) So if I hold down `command` on a Mac keyboard or `control` on Linux or a PC, you can see that these are actual links.

[11:43](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=703) So if I want to find out more about what this is, it will take me to an nf-core webpage about this specific type of test.

[11:57](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=717) We’ve updated this, so the soft-link didn’t work, but you can see long documentation about all of the different nf-core tests that we have here.

[12:06](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=726) This is something we’ve just rebuilt extensively in this release. That’s how all this automatic linting works, but you can go to the nf-core website to figure out all of the reasoning behind each test and how to fix it.

[12:17](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=737) Some of them are not super obvious, so we can go through it step-by-step.

[12:24](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=744) There’s some stuff here that is brand new for this release, which is nice to talk about as well.

[12:29](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=749) The bottom-up we have some stuff. We have a lot of tests passed, which is great.

[12:35](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=755) By default, there’s a flag if you actually want to see every single test pass. But most of the time you don’t care.

[12:58](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=778) When you create a new pipeline, you will see lots of things here. So these are not going to give you a little red cross but it’s good to cut that list down.

[13:12](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=793) In this version of the tools, we have the ability to ignore link tests, which is a very much requested feature.

[13:20](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=800) So I actually have a file in the root of my pipeline here. If I do `ls -a` to reveal all the hidden files, and then I do `cat.nf-core lint`, and you can see that I’ve specified in this config file and tools that the lint test (called files unchanged), should ignore this specific file.

That’s great because I’ve edited this and it would be failing anyway. But instead, it’s recognised this and ignored the file from this test.

[13:45](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=825) You can fully ignore any lint tests now and you can customise which parts of certain tests such as which files to ignore on that level.

This is really powerful, especially if you are using nf-core/tools but are doing your own pipeline that has got nothing to do with nf-core.

[14:02](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=842) You can see the problem here is that something failed.

[14:06](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=846) There’s a test called `files_unchanged`, which checks whether files match the template, and it says that something has been edited in this file.

[14:15](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=855) So basically this file shouldn’t be touched, so that’s wrong.

[14:19](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=859) One of my favourite new features for this version of the release is that `nf-core lint` now has a --fix, which is going to save everyone so much time.

So before, you would have had to figure out, read the documentation, figure out, go look in the template, see what it should look like and then try to fix it yourself.

[14:34](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=874) So now I can just say that I’d like to fix this test called `files_unchanged`. I can run linting again, and this time it’s going to not only spot them but also fix them.

[14:47](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=887) I don’t have any tests failing now and you can see that I fixed this test. If I do get status, you can see it’s actually modified one of these files and if I could get this, you can see that I had this extra line that I had written in here.

[14:56](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=899), and it’s deleted.

[15:01](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=901) So now it’s matching the template and everything works.

[15:06](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=906) So now if I commit this, this little red cross will hopefully turn into a green tick for the nf-core/tools.

[15:18](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=918) You can see that I could have also fixed the `Conda` updates, which is just a nice and fast way to update all of the different `Conda` packages, but I chose not to do that this time.

[15:29](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=929) This only works when you have no changes on your `git` repository and that means it will make whatever changes it can.

It can be quite aggressive, but if you don’t like it, you can always undo it by checking out the old code.

[15:44](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=944) So if I do update this test, it gets rid of all those warnings about Conda packages.

[15:56](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=956) You can see it has modified my environment file and updated all these packages.

[16:05](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=965) I can always do `git checkout environment`, and I haven’t lost anything.

[16:11](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=971) OK, so we go back here. Let’s see if this latest test works.

[16:27](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=987) Fantastic, we’ve got the little green tick across the board! That’s what we like!

[16:32](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=992) There are some little extra things that might be useful. You can see there’s a button up here [16:34](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=996) for some of these tests. There isn’t for the `markdown linting` though.

[16:41](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1001) It has automatically saved the verbose log file from the linting run, so if you can’t figure something out, or you need extra information, try downloading this [16:49](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1009). It’s a massive log file that’s spitting out debug messages about how nf-core/tools is running. That might help you debug a little bit.

[17:00](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1020) Anyway, we’ve got all green tests there, so let’s open up a pull request.

[17:04](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1024) We go to the head pipeline, go to my fork, create pull requests.

[17:19](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1039) Don’t just delete all this, you should actually fill it in, but you know this is a live demo and I’m going to create a new pipeline.

[17:25](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1045) You pull request.

[17:29](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1049) So what I’ve deliberately put in here was a merged conflict because it’s something that happens quite often.

[17:35](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1055) I had said that tests can run on different event types, so the ones we’ve been looking at are push pull requests, and there’s also a pull request.

So when you open a pull request, you will have tests running on both push and pull requests, but if you have a merge conflict, nothing will run on the pull request. So you need to fix all the merge conflicts before the test will run.

[17:59](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1079) So while I was talking, some magic happened and we had github actions’ automated comment pop-up.

[18:06](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1086) I made another deliberate mistake when I opened this pull request. It is against the `master` branch. Now with nf-core/pipelines, we have only the stable released code on `master branch` because that’s the one that’s pulled by default.

[18:17](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1097) So it should have gone to the `development branch`, but it’s really easy to mess up because it’s a default.

[18:23](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1103) It happens a lot, so we get a comment that says, “Well hang on, one of the tests here failed.”

[18:31](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1111) There’s an automated test to check it, so I got a red cross and I got a comment saying this is what has happened and is what’s wrong. You don’t need to close and open your pull request again, you can just hit edit, change that to dev and we’re all good.

[18:46](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1126) So this comment won’t go away, but we can hide it and say it’s been resolved.

[18:52](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1132) This little red cross will also not disappear until I push a new commit. But I need to resolve the conflicts anyway, so let’s just do that quickly.

[19:01](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1141) If I go to the repository (and look at my remotes), I’ve got my fork set up as a remote called origin. I’ve got the main nf-core one set up as upstream.

[19:13](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1153) So I can do a `git pull upstream master`, and it tells me that I have a merge conflict.

[19:20](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1160) So if I hop into vscode, which is here, you can see the `merge conflicts` and `CHANGELOG`.

Sure enough, there are two lines, which have been added kind of in parallel, and git does not know how to merge in both. But, it’s just a `CHANGELOG`. We want both of these [19:33](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1173).

[19:34](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1174) So I can click accept both.

[19:37](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1177) There are some differences in the markdown styles here, so hopefully if I hit save, it will solve that.

[19:43](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1183) I can go back and `fix merge conflict`. This will do several things. It’s fixing a merge conflict so the pull request test will run. It’s also pushing a new change after I changed the target branch to dev.

[20:00](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1200) So the target branch test should now be passing, and we’ll see how much else is green.

[20:08](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1208) This is the point where I suspect we won’t get green across the board because of running on the development version of nf-core. I think we didn’t manage to scrape up every single test, but hopefully by the time you run this you will get all green tests at this point.

[20:39](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1239) So what I haven’t talked about so far is running workflow tests.

[20:43](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1243) So every nf-core pipeline has a way you can run it where you do `nextflow run pipeline name` and then you do `-profile test`. That runs the pipeline with a very tiny test dataset, which is downloaded from the web dynamically by an expert.

[20:57](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1257) So now basically, every single time you push or open a pull request, this test profile will be run.

[21:05](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1265) That just run the Nextflow pipeline and checks that it doesn’t crash and picks up a bunch of kind of potential problems here.

I’ve forgotten here that these tests are also quite intelligent. They check whether you have changed anything in the environment (the Docker file or Conda). If you have, then it will build a new Docker image before running the pipeline. If you haven’t, it will just pull it from Docker, which is a lot quicker.

[21:31](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1291) So we have to wait for this now to build a new Docker image, which might take too much time. But after that, you can see it’s going to install Nextflow and will then run the pipeline and check that it exits with a successful exit code.

[21:42](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1302) So that’s something you need to make sure is passing, and that can be run locally too.

[21:51](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1311) OK, so you see all your lines green and when this is all done, it will hopefully all be green as well.

[22:02](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1322) This is also something I’ve fixed for this release. So what you should start getting now is that when nf-core runs on pull requests, you should automatically get a comment.

[22:15](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1335) It would be very easy to not go into a log and see all the warnings, but we can get a summary here and see how actually we had some warnings. So this helps make those more visible and hopefully you can kind of solve those warnings.

[22:28](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1348) You’ve got some nice rendering of those same results.

[22:36](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1356) So when that all goes green, you get this screen across the board.

[22:38](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1358) You’ll need to get a pull request approval and then you will be able to take this and merge it because everything will be green.

[22:44](https://youtu.be/U9LG_mMQFMY?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1364) So that’s your continuous integration tests. There’s a Slack channel for linting, and that’s where I’d go to get help if I’ve encountered any problems with a pipeline and nf-core lint tests. Don’t suffer in silence.

</details>
