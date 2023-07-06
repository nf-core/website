---
title: 'Bytesize 16: Modules test data'
subtitle: Kevin Menden - QBiC Tübingen, Germany
type: talk
start_date: '2021-06-08'
start_time: '13:00 CEST'
end_date: '2021-06-08'
end_time: '13:30 CEST'
youtube_embed: https://youtu.be/QXfAerydAT0
location_url:
  - https://youtu.be/QXfAerydAT0
  - https://www.bilibili.com/video/BV14y4y1u7ns
  - https://doi.org/10.6084/m9.figshare.14748438.v3
---

# nf-core/bytesize

Join us for an episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 16: Modules test data

This week, Kevin Menden ([@KevinMenden](http://github.com/KevinMenden/)) will present: _**Modules test data**_

This will cover:

- How to use existing and add new test data for nf-core/modules

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:53](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=53) Thanks for the introduction. Today we will talk about nf-core/modules and test data for these. Specifically, I will be talking about how to use the test data that we have, to create tests for your modules. I will also cover what to do if you don’t find the data you need for the specific module that you want to add to the nf-core/modules repository.

[1:24](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=83) So I would like to start with a quick recap about tests for nf-core modules that have been covered in previous bytesize talks ([Bytesize#5](https://nf-co.re/events/2021/bytesize-5-dsl2-module-development),[Bytesize#6](https://nf-co.re/events/2021/bytesize-6-dsl2-using-modules-pipelines)).

[1:26](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=86) ..and there’s also a lot of documentation online. So be sure to check those out. But as a reminder, for every module we have on the nf-core/modules repository, we have unit tests for Docker, Singularity, and Conda. These tests consist of several little scripts. The main script for the test code is located in `modules/tests/software/fastqc/main.nf`. There is a separate directory for every module, and there is a `main.nf` file that contains the main code. There’s also a `test.yml` file, which we will look at in a second. Then there is the `pytest_software.yml` in the tests config directory (`modules/tests/config/pytest_software.yml`). This file contains short entries for every module, so you need to ensure that your tests are in this file so that the tests are executed by GitHub actions.

[2:36](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=157) Let’s now have a look at one of those main.nf files. So this is one from the `fastqc` module as an example. As you see, it’s a simple file, you include the module `FASTQC` from the software directory, and then define a quick workflow where you basically do two things, the first of which is to define some input data, which is what we will be talking about today. So you see here, you define a file, then link to the test data parameter. Then there are a couple of keys to define the test file that we will look at in a bit. Then you just run the module on this input data and get an output. That’s the test.

[3:25](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=205) The next thing that comes with this is the `test.yml` file, where you define some information about the test. In this case, the name of the test is `fastqc single-end` i.e. fastqc for single-end data. Then you define the commands for the test, and a tag for the test, which is used to run the rest by GitHub actions. Then also importantly, we have the output files defined. In this case, the module produces both a .html file and a .zip file as outputs, and usually there are also md5 sums in here so that we can see that the output of the module is the same. That’s not the case here, but it doesn’t matter.

[4:14](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=254) So those were the files that made up the tests. Now, you might wonder how you actually know these dictionary keys I showed you earlier. To define the file that we’d like to use for testing, we have this command file which is a Nextflow command to define a file. Then we have the test data parameter, which is basically a dictionary. Within that, we look for the `test_1_fastqc_gz` entry using the `sarscov2` and `illumina` keys. If you’re wondering where to get these keys from..

[4:54](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=294).. For this the most important file is the `test_data.config` file, which is stored in `modules/tests/config/test_data.config directory`. So this file is used to define all these different test data files that are available. It has a link to the test data directory, which is the nf-core test dataset repository, and specifically the modules branch of this repository where we store all the test data. There are links to specific files in this directory, and these links are again linked in this dictionary. This is just a small example of this whole test data parameter dictionary. The reason we do it like this is so that we are able to change these links at some point without changing the code in the `test_main.nf` files of the modules because you will only have these keys here and these entries are defined in your script, so you don’t really need to care about what’s happening on this side of the code. So this is where I can look for the keys and the data I need in the file.

[6:12](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=372) So where can we get this data from? For this you need to look within the `modules` branch of the `nf-core/test-datasets` repository (<https://github.com/nf-core/test-datasets/tree/modules>). This repository is where we store all the datasets that we use for testing. There’s a `README` file, where we try to explain these different test datasets. Please make sure to take a look at that. Then also go through the repository and go through the datasets you need, to ensure that they are the correct ones.

[6:52](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=412) We currently have mainly two genomes available, one is the sarscov2 and the other is a subset of the human chromosome 22. The reason for this is that sarscov2 is really small, so we can really use the whole genome, and use all kinds of files for this genome without having large datasets. This is important for two reasons, first, we want to have small datasets to test, so that the test repository does not explode and the tests are also smaller for small datasets. For the human genome subset of chromosome 22, we again have a small dataset. Currently we only have genomics data available on the test dataset repository, and we also have some non-standard genomics data for specific tools that have been added by users.

[7:54](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=474) So how do you use that data when you write a test? The first thing you need to do is to identify the file you would need. To do that, as I previously mentioned, you need to visit the test dataset repository (go through the README) and take a look if the file is available, double-check where the file is located, and finally use that information to look for corresponding keys in the `test_data.config` file. You can use those keys to link the file in your test code. At that point, you should really be able to just run your test and finish your test code. That’s it then!

[8:32](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=512) Unless of course you can’t find the data that’s necessary for your modules. So of course we don’t have all kinds of datasets online and there are all kinds of bioinformatics tools available. Some of them may need specific data to run the modules. So it’s quite possible that you will find that the data you want is not available. This isn’t a problem, you can just upload the data yourself (this will also be helpful for us!).

[8:57](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=537) So first, it is important to read the instructions (in the `README` in the test dataset repository). We have a few lines about how to add new test data. Make sure to check this out and then think about a suitable minimal test dataset. The optimal thing here would be if you could find something that went with sarscov2 because that’s currently our favourite genome because it is small and there’s a lot of data available. So hopefully you can use this because then everything stays small and we don’t need to add a new genome. Then once you identify something, it is best to ask for opinions to check if it's suitable for the test dataset repository in the `#modules` channel on Slack. Once you’ve established that it’s OK, you can create a pull request to the `modules` branch of the test dataset repository and ask for reviews `#modules` so that you can eventually have your pull request merged. Once that’s done and your data is available, you need to add an entry to the `test_data.config` file in the `modules` repository. That’s again a separate pull request to have an entry merge, but it shouldn’t be a problem.

[10:18](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=618) Once that’s done, you can use your data just like any other test data that is available.

[10:26](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=626) So now for some guidelines for new data because we want to keep all the data in this test dataset repository somewhat similar and small. First, the data should be as small as possible, but of course big enough so that there is some output, yet the output does not need to be meaningful. Then, try and use existing genomes when possible, and if you add new data, try and adapt them to these genomes. So for example, if you would like to add human data, you need to ensure that it is compatible with the region of chromosome 22 that’s already there. Finally, if it isn’t possible to use the two existing genomes, we can have a discussion on Slack and add a new genome if necessary. We also like standardised names on the test dataset repository, the goal is to have the files similar across the different genomes, so that we could theoretically just swap the genomes and the test names would stay the same.

[12:06](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=726) Adding test data can be a bit of work because it involves first looking for the data, or needing to produce it, cut out certain regions that fit the human genome, and verify that it is the human genome. So it can be a bit time-consuming. It is however worth it because it would not only help get your module on the nf-core/modules repository, but it would also be useful for other modules in the future. That also makes it really useful for the community as a whole. We envision that this might also be beneficial for other communities as well since we hope to merge different test dataset repositories with those from other initiatives.

[13:25](https://youtu.be/QXfAerydAT0?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=805) So that was it from me. If you have any questions, please feel free to reach out to me on [`#bytesize`](https://nfcore.slack.com/channels/bytesize).

</details>
