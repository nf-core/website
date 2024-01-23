---
title: 'Bytesize 18: Development environments & workflows II'
subtitle: Maxime Garcia - SciLifeLab / Karolinska Institutet, Sweden
type: talk
startDate: '2021-06-22'
startTime: '13:00+02:00'
endDate: '2021-06-22'
endTime: '13:30+02:00'
youtubeEmbed: https://youtu.be/OF55x-FT5WE
locationURL:
  - https://youtu.be/OF55x-FT5WE
  - https://www.bilibili.com/video/BV1Wh411679X
---

# nf-core/bytesize

Join us for an episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 18: Development environments & workflows II

This week, Maxime U. Garcia ([@maxulysse](http://github.com/maxulysse/)) will present: _**Development environments & workflows II**_

This will cover:

- Git + GitHub configuration
- VSCode configuration
- Pipeline testing

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:34](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=34) Hi everyone, the talk today will build on [Bytesize#11](https://nf-co.re/events/2021/bytesize-11-dev-envs-workflows).

[0:45](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=45) I will go over my own best practices for working: the configurations that I use for Git and GitHub. I’ll also cover a little of VSCode Configurations and some local testing as well. This has also been covered by Phil in [Bytesize#11](https://nf-co.re/events/2021/bytesize-11-dev-envs-workflows), and my workflows aren’t that different from his.

[1:23](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=83) But I'd like to issue a disclaimer; what I’m about to describe is what works for me. It’s probably not the best, but it might be useful for some.

[1:37](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=97) So here are some of my best practices for working. I spend 20 minutes each morning on email and Slack; I would encourage you to use an email filter and for Slack I recently discovered that you can set your preferences to see all your unread messages at once. I really like that and it helps me save time. I also try to stand up for five minutes every hour and take a short walk around the apartment, drink some water, and stretch a little. I think it’s important. When I really want to focus on my coding, I avoid being distracted by muting all my conversations on Slack. It’s also important to take a break when you need one.

[4:52](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=292) To follow what was presented in [Bytesize#11](https://nf-co.re/events/2021/bytesize-11-dev-envs-workflows), I’ll cover how to set up your own fork on GitHub, configure your local clone, work on other people’s code and pull in updates.

[5:09](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=309) I use the GitHub command line. Let’s say that I’ve forked the pipeline. What I usually do is to append the name of the original directory to my fork, so that I can trace it back if I need to. If I’d like to check the code, I use the GitHub command line interface. I copy the path and go to my workspace folder and then clone it. An advantage of using the GitHub command line interface is that the remote is already configured. This makes it easy.

[6:26](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=386) Then if I’m doing some actual work. I use `zsh` on my terminal so that I know which branch I am on (I’m on `dev` here). When I start working on a new feature, I do `git checkout -b` on the branch I am working on and then specify my new feature. This takes me to that branch, and I work on my branch where I’m now going to modify a file. Let’s head to VSCode, take a look at the pipeline, then head to the CHANGELOG.md, and make some small changes. I’m happy with the modification. So now I go back to the command line and say `git add` followed by the name of the file. Then I commit with the message `git commit -m` and add the feature (“feat:”) to update the file name (`CHANGELOG.md` in this case). I specify “feat:” if I’ve done something, “fix:” if I’ve fixed something, “chores:” if I had to do something. Then if I want to push, I use `git push -u origin` with the name of the branch. So that’s it. My branch should now be available on GitHub. You can also automatically create your PR in the GitHub command line interface with `gh pr create`. That’s quite fun and it’s really useful to make a PR like I needed to do when I updated the social media images for the different pipelines. You need to specify the base branch `--base dev`, a title `--title` for example “test”, and you can specify a `--body` for the message for example “testing something, do not merge”. When I do that, it asks me where I’d like to make the repository, and I add that information.

[10:40](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=640) So if I want to make some tiny modifications, then I do that directly on the GitHub interface. Let’s look at updating an event on the nf-core website; specifically to remove a zoom link for one of the bytesize talks. So I just remove the link, and write a meaningful commit message. GitHub automatically names the branch with your GitHub login ID and patch - x where x is a number. I usually change that to reflect the branch I am creating a patch for, so in this case my patch would be called “maxulysse-bytesize18”. I then click on 'propose change', ask for it to be reviewed and create a pull request. So you can do a lot directly on the GitHub interface and it’s quite convenient.

[12:51](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=771) Now let’s see what I do if I want to keep the code up-to-date locally. I first do `git status` to check for uncommitted changes, and then `git pull` to check if there are updates. I also like to do `git fetch upstream` and `git merge upstream/dev`. If I make some changes and then do `git push`, then the branch will be up to date. If I now go to my branch (on the GitHub interface), then I see that it isn’t updated. So I click on `fetch upstream` and then `fetch and merge`. You can also do that on PR, and I would really recommend that there.

[14:56](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=897) Now for some best practice and useful plugins on VSCode.

[15:06](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=906) I organise VSCode is by having a big workspace folder, and if I want to work on several different things in multiple repositories, I have it all in the same workspace. But if I want to work with just one specific folder, then I open that workspace. For example the Sarek workspace which just has the code for Sarek. This allows me to focus.

[16:27](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=987) So that covers how I work with VSCode. I try to do things one project at a time. But you should do what works best for you.

[16:52](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1012) I don’t really have any specific recommendations for a plugin.

[17:40](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1060) So now let’s talk about how to run tests locally. I essentially do more or less the same as was described in [Bytesize#11](https://nf-co.re/events/2021/bytesize-11-dev-envs-workflows). Let’s use the sarek pipeline here. I run the tests that I’m doing in the GitHub command line, for instance, I am currently working on annotation on the DSL2 branch, so I run a test that’s specific for annotation. So I run the command line to see if it is working, then I wait for Nextflow to start. We will see this test running in the CI on GitHub actions. After last week’s talk ([Bytesize#17](https://nf-co.re/events/2021/bytesize-17-pytest-workflow)), I started working on adding tests for specific sub-workflow or for the whole pipeline. So if I look here at my test, I can run the same test separately with a specific tag. So in the command line, I go `pytest --tag annotation --kwd` if I want to run all the tests or `pytest --tag vep --kwd` for a vep test. It’s easy to have this set up and running once you have pytest installed. I also always run it in the command line first to check that it runs like I expect it to and then run it on pytest with the specific tags.

[21:12](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1272) I tend to listen to music if I am having difficulty focussing. We have a [dedicated playlist on Spotify](https://open.spotify.com/playlist/6LyhtB3bllSwNbK9iDNVgH?si=da2dd5f0b6134b7b) if you’d like to listen there. Also please feel free to add more music there if you’d like. It’s a collaborative playlist!

[22:05](https://youtu.be/OF55x-FT5WE?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1325) I would like to thank you for listening, the institutes and funding agencies that I am supported by, and all the institutions that contribute to nf-core. For more of a background on what I discussed today, please refer to [Bytesize#4](https://nf-co.re/events/2021/bytesize-4-github-contribution-basics), [Bytesize#11](https://nf-co.re/events/2021/bytesize-11-dev-envs-workflows), and [Bytesize#17](https://nf-co.re/events/2021/bytesize-17-pytest-workflow).

</details>
