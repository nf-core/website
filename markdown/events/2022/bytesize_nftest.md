---
title: 'Bytesize: nf-test'
subtitle: Edmund Miller, University of Texas at Dallas
type: talk
start_date: '2022-12-06'
start_time: '14:00 CET'
end_date: '2022-12-06'
end_time: '14:30 CET'
location_url:
  - https://doi.org/10.6084/m9.figshare.21695195.v1
---

# YouTube video

<!-- markdownlint-disable -->
<iframe width="560" height="315" src="https://www.youtube.com/embed/K9B7JRkMpQ4" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
<!-- markdownlint-restore -->

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: nf-test

This week Edmund Miller ([@Emiller88](https://github.com/Emiller88)) will share with us his impressions about nf-test from a user perspective. nf-test is a simple test framework for Nextflow pipelines.

<details markdown="1"><summary>Video transcription</summary>
**Note: The content has been edited for reader-friendliness**

[0:01](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=1)
(host) Hello everyone, welcome to today's bytesize. I'm Franziska Bonath, and I'm the host for today. With us is Edmund Miller, and he is going to give us his impressions of working with nf-test.

[0:16](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=16)
Okay. Good morning, everybody. Before we start, I did not write nf-test, I've just been using it. This was actually written by Lucas Forer and Sebastian Schoenherr, I'm going to butcher that from iMed. Okay, so what is nf-test? It is basically the ability to write unit tests for Nextflow workflows, which I looked at Phil's talk, is the second most requested feature in this year's Nextflow and nf-core community survey. So you can think of it as Pytest, but for Nextflow, if you're familiar with Python or unit tests, basically just a way to write either full pipeline tests or sub-workflow function tests and module tests, and we'll get into some of those.

[1:14](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=74)
So what that looks like. If you haven't memorized the Hello Nextflow script, here's a quick reminder. We start saying all of these different greetings in different languages, and we pipe them through, and then we view it lastly. And so what that looks like on the nf-test side is then we have a name of the test, this boiler plate up here, and then we have a script, and we're actually calling, I believe, the remote script here, and then a name of the test on this line, and then we have what we expect, such as the success of it, how many tasks we're going to have, which I think is also important as you're changing things, and then we have several asserts of what's happening in the standard out of these greetings, for example. And that's the power of what you can do with nf-test here, and I'll let your imagination run wild with that of all the things you can do with it from there.

[2:12](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=132)
So a little bit of background as to why we want this, and if you're not sold yet on that. So far, our testing practices until now have been testing our pipelines with GitHub Actions, and then testing our modules with PyTest workflow that's running on GitHub Actions for those. So with the pipelines and GitHub Actions, I've realized I was talking to Sateesh as we were starting to play around with nf-test that CI is meant to run your test framework, not be your test framework. And I think that's rung true of all that it's doing really is it's just checking the workflow ran without failing. And that's a huge step up from not checking anything, but then that starts to create several problems as well when we're doing that, such as a number of the tests, whenever that starts to go up, it starts to become more complicated to maintain. You start coming up with all these fancy matrices, things start breaking. It's hard to keep track of what's going on in the jobs per se. And you can see Sarek as an example of that and all of the various tests that they had and all of these and trying to maintain their huge pipeline.

[3:32](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=212)
To run the test locally then and repeating those as well in a fashion, you had to really translate from that workflow YAML and then convert it to a Nextflow command and then run those. If you wanted to run multiple tests, you had to do that a couple of times. This had a heavy reliance on CI runners as well for us to produce those rather than us using the CI test or the CI runners as a check and a settlement of an argument of whose computer it runs on correctly.

[4:09](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=249)
In our modules, we've been using Pytest workflow. If you're familiar with that, I'm sure if you've contributed, you've fought with Pytest workflow at some point. Running your tests should be as easy as running the pipeline. Pytest workflow, however, was not. And so we created `nf-core modules test`. And that's a great utility. And it's been a great way for people to test things locally and rerun them. But it's a lot to maintain on our side and change and update it and it keeps us from wanting to progress forward on those. We want to keep the infrastructure the same so that we can focus on actually writing the Nextflow modules.

[4:54](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=294)
There's a couple of problems though with Pytest workflow. As great as it's been, I've loved it. And it's been a huge step up from nothing and being able to actually specify tests and lay them out. The issue is that we have no native support for Nextflow profiles, which is a real power of Nextflow. So we can't multiplex Docker Singularity Conda. And yes, we've come up with several hacks to get around that in CI and locally. But it's really difficult to reproduce locally and explain to new users of like, okay, well, Nextflow lets you do `-profile`, but then here you have to specify it before the workflow. And it's a lot of caveats for a new user and a new contributor. Whereas nf-test allows you to just do `-profile` just like Nextflow would. And so it's a lot easier for new users to say, okay, I want to test with Docker and you jump right off to that. There's not as many gotchas.

[5:51](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=351)
The other issue that I really had with Pytest workflow is the pipeline output was tough to find locally. You had to go dig through "temp". The directory changed every time. You can really go figure out where your output was for your Nextflow workflow. We wanted to convert Pytest workflow or convert all the CI to Pytest workflows. Why were we putting that off across the board? We had done it in a few places and Sarek had done it in all of it's tests as well and done lots of extensive work. But it was a lot of manual work to get all the md5sums for the files that we wanted. And the testing stopped there if you can test for all the files and you can test for md5sums and you can do contains and other stuff as well with Pytest workflow. But it was a lot of manual work to get all those there.

[6:49](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=409)
We can build new infrastructure like we did for the Nextflow or the nf-core modules and start to build all that and box ourselves back into a corner as well for this rather than just relying on a different framework to do that. The other big complaint that I had for Pytest workflow for pipelines was I could never get local modules working with the bin properly for the mocks that I wrote. With nf-test, it just works out of the box.

[7:23](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=443)
After all of that, now a quick intro to nf-test. All you need to do to get started is a simple curl, which if you've installed Nextflow at this point, you're pretty familiar with. It's also Java and groovy based and runs with Nextflow. We're hoping to have the nf-test folks back on for a more in-depth technical talk on this. We just want to get started on it. So the next command that you'll run is `nf-test init`. This just sets up the testing config for you and creates a testing directory and a default test config. You can play around with those. We haven't come up with the perfect best practice yet for what we want for the nf-core standard to be on that.

[8:12](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=492)
The next thing that I think is really, really powerful is they already have a `generate` command. So you can generate all of these tests for your pipeline and workflows processes. And also really importantly, we haven't tested the full ability of it's functions. Lastly I think that'll really start to change things, especially as we were starting to change to checking our sample sheet with groovy, that might help us specify those and lock some of that down. Lastly, you just have to run nf-test after writing all of those tests. That will generate all of the boilerplate that you need for basic smoke tests of all the processes, workflows, functions, et cetera.

[9:02](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=542)
Now a quick little walkthrough of some code. This is just a quick test that I wrote in demultiplex, and I don't think the highlighting's working yet. Oh no, there we go. What's really cool here is, we can start to specify params. So if you've gone through and done this in your workflows and thrown it into GitHub Actions and thrown this into a matrix, you know, you start to feel limited after you get to one or two, or you start to really specify out an entire test and it's just not very familiar. I believe you can parameterize these as well. So you can change your inputs and change various params for your workflow.

[9:50](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=590)
This is just some boilerplate. It boxes each workflow into its own separate execution directory, I believe. Phil and I have had some interesting findings with that. Again, you can check the size. Now what is the really powerful part here is this snapshot feature. What this does is, it basically pulls in the md5sums for you of all of these files. So we have things like metrics, the run manifest here. We have our fastqs that are zipped up and you can specify all these different paths. This actually isn't a perfect example, because you can just list a directory, I found that really powerful. You can just point it there and then it'll go through all the files. What that looks like over here is that you get this snapshot file.

[10:40](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=640)
What this does is it just creates a quick little snapshot with the md5sums of all those files that you specified. And you can do standard out, standard error, et cetera as well on those based on the workflow. But we found that to be a little hairy for right now. As you can see, this is automated. And whenever you rerun your test, it checks for these md5sums and will tell you. The other thing is it's really easy when you know you changed some stuff or you know these are going to change. You just run it with the `--update-snapshot` flag and it updates all the snapshots for you automatically. Then you just commit this file and it will check again against them.

[11:21](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=681)
If we can go back. For files that you can't md5sum, you can also assert these files and check that they exist as well. That's one thing. There's lots of other powerful features as well. You can do regexes and you can also start to build plugins for it, such as checking for FASTA files. I think there's plans for a VCF checker as well. So I think that'll be really cool to see what comes out of that.

[11:58](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=718)
Let's talk about some good testing practices as we start to adopt this and start to adopt more testing. I know we have great testing practices in the nf-core tools. So there's things like test driven development and that's where you write your test first and then you write your code after. There's also acceptance driven test driven development. That's where you start to write your acceptance criteria first and then you write your tests and do your development. There's behavior driven development, where you start to talk to the team first and then you go through those and start to write your tests and then do your development after that. There's also pair programming, which is where you write the code, someone else writes the tests, you come together and work on it together. Lastly... just kidding, we're not going to talk about any of those today.

[12:57](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=777)
Those are all great and lots of great things that you should go read about. But I think as a community, we have more important things to focus on. What are some realistic testing priorities? We should probably first convert our CI tests to pipeline tests. I think those are great to reproduce locally and quickly. Next, we can then add more pipeline tests for various params and pathways that you might want to make sure don't break, I think most people are just testing that all their aligners work and then hope that everything downstream works, whereas you might want to check all of your skips, for example, and make sure that you get the right amount of processes and that certain things aren't getting run.

[13:44](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=824)
Another great practice that I hope that we can adopt is to add tests when fixing bugs in PRs that prove that this bug is fixed and we're not going to come back to this bug and checking for it in the future so that we don't end up recreating it. Lastly, I think testing your local modules can go a long way. I think we have a lot of custom scripts that we could do some tests on and make sure that we're getting the expected output out of those to verify our scientific results. So the TLDR is let's avoid regressions and move forward with our tests.

[14:22](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=862)
A couple of us have been toying around a lot how a roll out plan could look like. This is just a quick update on this and is no way indicative of exactly what's going to happen. First, we want to start with pipelines as a proving ground because I think that's where we had the most ground to gain compared to being stuck with just GitHub actions. So it's already in methyl-seq and demultiplex. If you're looking for some quick examples and you want to implement this today and switch those over, just that step one, we're going to need to update the template. I've got a PR started for that, but we need to settle on some best practices or how we want to do things so that as you switch between nf-core pipelines and contributing, things are familiar.

[15:07](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=907)
We need to do some module infrastructure prep and figure out how we can convert those and run all these tests and what that's going to look like. We're waiting for them to come out with tags for nf-test. And then we'll update the modules as they get changes, hopefully, on those. That will be a roll out as you make a change to a module. You'll update the tests from the Pytest workflow to nf-test followed by a final push to convert all the modules that didn't get any love at a hackathon. And with that, I'll take any questions.

[15:43](https://www.youtube.com/watch?v=K9B7JRkMpQ4&t=943)
(host) Thank you very much. I am now allowing people to unmute themselves, also to start their videos. Really nice talk. Are there any questions in the audience? I don't see anything.

(question) Is there any logo for nf-test?

(answer) I don't think so. You should make an issue for it.

(host) No logo. That can't stand.

(question cont.) Yes, but definitely logo, stickers. Yes, Marcel, you're right. I want stickers as well.

(host) Okay, but if there are no other questions, then thank you again, Edmund. And I would also like to thank the Chan Zuckerberg Initiative, as usual, for funding these talks. If you have any questions, you can always go to Slack, ask your questions there. I guess for this one, it would be in a general help channel. If you have questions about nf-test or about any testing question. So thank you again.

</details>
