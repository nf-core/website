---
title: 'Bytesize 17: Pytest workflow'
subtitle: Edmund Miller - University of Texas at Dallas
type: talk
startDate: '2021-06-15'
startTime: '13:00+02:00'
endDate: '2021-06-15'
endTime: '13:30+02:00'
youtubeEmbed: https://youtu.be/pjhscKyWH74
locationURL:
  - https://youtu.be/pjhscKyWH74
  - https://www.bilibili.com/video/BV1k64y1r7WR
  - https://doi.org/10.6084/m9.figshare.14784855.v1
---

This week, Edmund Miller ([@edmundmiller](http://github.com/edmundmiller/)) will present: _**Pytest workflow**_

This will cover:

- What [`pytest-workflow`](https://pytest-workflow.readthedocs.io/) is and what we use it for
- How we test nf-core/modules
- How we test nf-core pipelines

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:52](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=52) Let’s get started with pytest workflows.

[1:00](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=60) We’re first going to cover pytest workflows and how we use it in nf-core/modules, then we’re going to duck behind the curtains of CI and understand how that works with the modules and how some of that is automated and generalised. Then we’re also going to talk about running pytest workflows on nf-core/modules locally as well, in case you have a lot of changes. Then finally, we’re going to get a little sneak peek at the testing pipelines with pytest workflow, and take a look at what’s coming in the future.

[1:32](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=92) So first an introduction to pytest workflows; what it is and why we’re using it. Writing workflows is really difficult and that’s why we are all here in a community called nf-core. Testing if they are correct is even harder; testing `bash` scripts or other code can be quite flawed. For instance, what are you going to use to test it? Is there a bug in the pipeline or is it in the test framework? So before we just used CI and ran the commands manually to check whether the pipeline even passes and doesn’t end with an error. Pytest workflows aim to make testing as simple as possible so you can focus on debugging your pipeline.

[2:20](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=140) So how are we using it in nf-core/modules? Here are some benefits of using a testing framework. We can version control and collaborate on our tests that we’re running on modules and pipelines. We can use cool things like `git bisect` to find out what commit broke something. This allows us to increase the reproducibility of our tests as well and it also allows us to run the test locally quickly without trying to copy commands from a CI file and convert them from github actions quirks to local quirks as well. It also allows new users to have a smoother learning curve in testing actual modules as well.

[3:04](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=184) So this is the quick anatomy of what a pytest workflow test looks like. If you’re not familiar with pytest, it is a very popular testing framework for the python language. Pytest workflow is a plug-in for that. It picks up special .yml files in the test directory. This is what those `.yml` files look like. So first you have the name of the test, and you can give it whatever name you’d like; something that would best describe it. Then you have the command that’s going to be run. In this case, for the example test, it’s going to run `touch test.file`. Then it’s going to check for the file in the path text.file after that command has run. The test will fail either if the file isn’t there or if the test doesn’t exit with a zero exit code. You can also specify what exit code it should have. So we can have things like tests it should fail because we’re testing for checks in that and various other things. I suggest that you look at the docs in detail.

[4:17](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=257) So let’s go through an example of `gunzip`. This is straight from the nf-core/modules, so you can reference that later. We have our name, our command which is `nextflow run` and then we’re changing to the directory that is for `gunzip`, and then we’re just going to use the entry point in the test file. There’s a `main.nf` in each test, and those have the Nextflow scripts that call the modules, run them based on that, and have some outputs on the test data. The talk that Kevin gave last week ([bytesize#16](https://nf-co.re/events/2021/bytesize-16-module-test-data)) refers to that. Then we’re just going to use the config in the `test config`. These tags here are also important so that we can tag the tests with different things, like `gunzip`. I’ll get into how you can use that in the CI later and also how you can use it locally. The other cool thing about this is that you can also `md5sum` the outputs and that is a cryptographic hash of the file. What it’s doing is that it’s checking for the integrity of the file each time it runs the rest. So then we can confirm that these are reproducible and we’re doing the same thing over and over.

[5:38](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=338) So luckily, you don’t have to know how to do any of that. You just need to know `nf-core modules create-test-yml`; this automates the creation of the .yml file and the md5sums as well. It also runs the tests for you for the first time, so all you need to do is write the test in Nextflow.

[6:10](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=370) This is an example of how Nextflow would do it. So first we put a test file in there; the STAR_ALIGN test workflow, which is just a main.nf in the test directory. You can refer to the nf-core/modules for this as well. So we just import those modules, have different parameters, we use the test data that Kevin talked about ([bytesize#16](https://nf-co.re/events/2021/bytesize-16-module-test-data)), and then just run the workflows in Nextflow.

[6:49](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=409) Then we’re just going to call `nf-core modules create-test-yml`. It will prompt you for a few things; you need to add some default values like `star/align`, define where you’d like the test.yml output etc. It suggests some defaults based on what you feed it as you see here in the parentheses. Hitting enter will probably be good. Then it’s going to look for the test workflow entry points as well; it picks those up automatically. It then creates a test name for you that follows our standardisation practices and then it’s going to come up with a test command for you - the test tags. You can also change any of these, add to them, or even go in and edit the .yml manually. It will also ask you for the test output folders and then run the test here. It will also pick up any other test that you have written; you don’t have to stick with just one, you can write multiple tests.

[8:06](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=486) So let’s peek behind the curtains of the CI.

[8:10](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=490) There are a couple of pieces to it. Luckily, GitHub actions has a beautiful way of working together and we can reuse other people’s code and actions that they’ve created. We use the dorny/paths-filter that checks for the changes that have occurred in the pull request and in the pushes, allowing us to only test things that have changed in the code. Based on these changes, we then create a matrix of the jobs, and the tests are run in the containers based on the matrix, and then in `pytest-workflow` and linting against the tags that we pick up in the first step. The logs get uploaded if something fails, so we can download and look at them to understand why they have.

[9:11](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=551) This is the part that is checking for changes in the CI. I’d like to draw your attention to the filters part, and specifically to the `tests/config/pytest_software.yml`. This is something you need to edit when adding new modules.

[9:36](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=576) You just go in there and you add in the tag that you want to use and I’ve shown the same one as the previous example. So what this does is that it asks it to pick up anything on this defined path and run the tests. If we make any changes to the module or any changes to the test, we need to re-run these tests.

[10:10](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=610) So now for the matrix... We support the current Nextflow version (we can support other Nextflow versions in the future too). It passes the tags for the first step based on the check for changes and then creates them based on a profile so we’re basically multiplying each of these by how many there are, and the tags one is variable based on how many tags there are. That’s created a problem for us in the past because it creates more than 256 jobs, GitHub actions kicks us off and doesn’t run them because the matrix is too large. I’ve tried some preliminary stuff to alleviate this and basically stop it before it gets to 256, but I haven’t yet got it to work. So if anyone has any ideas, that would be helpful.

[11:18](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=678) So then we’re going to run the tests. This is how simple the end of the CI is. We need to install python, pytest workflow, the profile, Nextflow, etc. All that it’s saying here is that it is the temp directory for singularity and we use the matrix tags here, so that it’s going to run back through each time and run it with different tags. Then we have `--symlink` and `--kwdof`. The latter deletes any working directory that hasn’t failed and cleans up space automatically.

[12:09](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=729) Any logs that fail here are uploaded as well. These are all just pytest workflow niceties that capture the standard out and the standard error for the workflow that’s running. We also upload the work directory (you may want to run that locally if you want to dig deeper since it might get a bit large otherwise). The logs are just named based on the matrix as well.

[12:43](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=763) We also lint the modules that use the same concept of checking for any changes, and run `nf-core modules lint` on the tags of those modules to save CI time.

[12:59](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=779) Now let’s cover how to run pytest workflows locally.

[13:05](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=785) It’s quite simple to get up and running. All you need to do is `conda install pytest-workflow` (you can refer to the pytest workflow docs if you’d like to do anything fancier). Then what you do is pick a PROFILE (this one is `docker`, but you can use `singularity` or `conda`) also depending on what you would like to run. Then you just call your tag of whatever module you are working on, and it will pick up all the tests for those. You can also use the name of the test as the tag; you just put quotes around it and it counts as a tag as well. Then we need `--symlink` in there for running it on the modules. You can also change the `--basetemp`, usually it puts it under `/tmp`, but maybe we want it under our scratch directory and not in that local space.

[14:07](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=847) So here’s a little sneak preview as to what’s coming in the future. We’re testing pipelines with pytest workflows.

[14:19](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=859) We can verify the expected outputs and their integrity, and this has been a longstanding issue. The CI can be generalised; it’s not so complicated and our logic isn’t in the CI. Our logic is in the test locally. This reduces our reliance on GitHub as a platform and allows us to be more decentralised if required. We can also only test sections of the pipeline that have changed. For example if we make a change in a subworkflow, we can just test that. This means that we save time on waiting for CI to run.

[15:11](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=911) Here are some example PRs for anyone who’s interested. If you’d like to implement this in your pipeline, get in touch via [Slack](https://nf-co.re/join).

[15:28](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=928) Here’s an end-to-end example that tests the entire application. So as you can see, this is a quick example from the `nf-core/rnaseq` pipeline (from running the default pipeline). We have the `nextflow run main.nf`. We can tag it with default and maybe we have some other default tests as well. We can also have the different files where we expect the outputs to be. This is where we benefit from the added complexity of pytest workflows. I picked examples that didn’t have md5hashes, but we can hash those in and verify the reproducibility. We can test the sub-workflow as well. We can input different options and test various options with your subworkflows. The advantage here is that you don’t have to run the entire pipeline and wait for Nextflow to spin up with these large datasets. Rather, you can just test your sub-workflows and work on those, and then plug it into the entirety of the pipeline when you want to test it. It makes for a very quick local development workflow.

[17:04](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1022) This is what it would look like for testing the sub-workflow. You see that it’s similar to how it’s running in the `nf-core/modules`, we have the various tags. We can use those to run the tags based off changes. As you can see, these are set to pick up anything in the `markduplicates` test in the sub-workflow file. If anything changes in the `nf-core/modules`, you will need to re-run those tests to ensure that the new version works. Of course we can have our expected outputs up here as well.

[17:56](https://youtu.be/pjhscKyWH74?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1076) That’s it. Get in touch if you have any questions.

</details>
