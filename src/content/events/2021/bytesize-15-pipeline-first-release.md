---
title: 'Bytesize 15: Pipeline first release'
subtitle: Alexander Peltzer - Boehringer Ingelheim Pharma GmbH & Co. KG, Germany
type: talk
startDate: '2021-06-01'
startTime: '13:00+02:00'
endDate: '2021-06-01'
endTime: '13:30+02:00'
youtubeEmbed: https://youtu.be/1OwkTd-P5pQ
locationURL:
  - https://youtu.be/1OwkTd-P5pQ
  - https://www.bilibili.com/video/BV185411M71d
  - https://doi.org/10.6084/m9.figshare.14717307.v1
---

This week, Alexander Peltzer ([@apeltzer](http://github.com/apeltzer/)) will present: _**Pipeline first release.**_

This will cover:

- Which steps are required to release a pipeline on nf-core
- How to make a first release of a pipeline

The talk will be presented on Zoom and live-streamed on YouTube.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:31](https://youtu.be/1OwkTd-P5pQ?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=31) Thank you for the introduction. I’ll be talking to you about pipeline first release, which will cover how one can release a pipeline that is built using `nf-core/tools`.

[1:06](https://youtu.be/1OwkTd-P5pQ?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=66) So what you need to do is to first prepare the pipeline prior to its first release. This is different from a continuous release, which is what is done on existing pipelines that are being further developed. The first release of a pipeline that has never been released before is handled differently.

[1:37](https://youtu.be/1OwkTd-P5pQ?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=97) So let’s jump in! We’ve covered certain prerequisites in previous `nf-core/bytesize` talks, and there are some useful links on the [slides](https://doi.org/10.6084/m9.figshare.14717307.v1). These range from how to sign up to become a member of the community (an obvious prerequisite for releasing a pipeline on nf-core) to how to create a pipeline from the nf-core tools template, add test data to the test data repository or the modules repository (test data for the DSL2 version pipelines), and add the pipeline to the nf-core organisation. So I’m now going to assume that you’re familiar with these steps, and are ready to release your pipeline.

[3:01](https://youtu.be/1OwkTd-P5pQ?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=181) So the very first step is to check that your pipeline follows the guidelines. This includes it passing all the continuous integration tests (unsure of how to do this? Check [bytesize#7](https://nf-co.re/events/2021/bytesize-7-nf-core-ci-tests)). You can go to your pipeline page, check Github actions, see if all those tests have passed, and then strike that off your list. (Make sure to check and see if nf-core linting has issued warnings; linting errors can result in CI tests failing). One way to ensure that everything passes is to have the most up-to-date version of containers such as Conda. It’s also essential to make the `master` or `main` branch the default branch. A link to how you can do this is also on the [slides](https://doi.org/10.6084/m9.figshare.14717307.v1) accompanying this presentation.

[4:13](https://youtu.be/1OwkTd-P5pQ?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=253) Now that that’s covered, let’s talk a little bit about nf-core/tools. There are tools within nf-core/tools that can help you prepare your pipeline for its very first release and subsequent releases. It is not very different from what you do if you have released a pipeline before. One thing you should do is to bump the version numbers on your `dev` branch, for example you can run `nf-core bump-versions 1.0.0`, for your very very first release. nf-core/tools will take care of bumping steps in your pipeline code. If you’re keen on learning how to do that, there’s also a [tutorial](https://nf-co.re/tools/#bumping-a-pipeline-version-number) that explains how nf-core bump versions work. The same applies to the `CHANGELOG.md`. Ideally, the `changelog` should list what you’ve been working on in your pipeline. For the first initial release, people tend to explain what the pipeline is capable of doing, describe the features present in version 1.0.0 etc. You can actually add all that to the `CHANGELOG.md` file. After this, you can open a `pull request`, from the `dev` branch to the `master` or main branch, and hope that it passes all the CI tests and criteria for the template requirements that we have for nf-core. Most pipelines in nf-core use [semantic versioning](https://semver.org), and a lot of the pipelines also use code names for pipeline releases - nf-core/eager for example uses the names of Swabian cities. You can also generate other code names, and this is entirely up to you whether you would like to use a code name for your pipeline at all. This link might be useful if you do: <https://www.codenamegenerator.com/>.

[6:37](https://youtu.be/1OwkTd-P5pQ?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=397) So then you can start with the review phase or the so-called fake pull request, which is a thorough core pipeline review by members of the nf-core [core team](https://nf-co.re/about) or more experienced pipeline developers. They will go through all the code in the pipeline to ensure that it adheres to the community guidelines. What we do is to first open a fake `pull request`, and then do a review against the very first `commit` in your pipeline code. This ensures that what’s been changed in the entire pipeline code, documentation etc can be checked. This is a bit more thorough than a simple `pull request` that you would use when updating the pipeline from version 1.0 to version 1.1. Please ask for help here because making a fake `pull request` can be complex, and it’s only rarely done. Ask for help on Slack `#request-review`, `#new-pipelines` or just ping the nf-core core team. Anyone in the community can actually take a look as well, since this involves the whole community. People can have a look and add their ideas or if they spot errors or potential issues that can be dealt with prior to the first release. It’s also up to you to fix all the bugs in the pipeline. Once the core-review team is satisfied with the changes and considers the pipeline to be within the scope of the nf-core guidelines, we will close that `pull request`. It doesn’t need to be merged since it’s a fake one, Its purpose is to allow initial review of the core pipeline.

[8:54](https://youtu.be/1OwkTd-P5pQ?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=533) Now after all that, you’re ready to make a release `pull request`! You should have had two reviews for your `dev` to `main` branch `pull request`, so you can actually merge that to the `master` branch, and after that it’s really up to you how you’d like to create your new release. There’s a really nice [guideline on how to do that on GitHub](https://docs.github.com/en/github/administering-a-repository/releasing-projects-on-github/managing-releases-in-a-repository) itself. It is essential that you use the same version as in the nf-core bump version. So if you bumped all the code versions to 1.0.0, then you should use 1.0.0 without a prefix such as ‘v’ for version, to avoid conflict issues with Docker images and Conda environments etc. It’s also a good idea to have a code name for the release description, and that is something you could copy from your CHANGELOG.

[10:20](https://youtu.be/1OwkTd-P5pQ?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=620) After that you’ve made the very first release of your pipeline! Now, there are some additional optional steps. You could for instance ask someone from the core team to get you a document object identifier (DOI) that allows easy citation. This is something that a lot of pipelines currently have, the nf-core/atacseq pipeline has a DOI badge in the README. There is always a to-do in the main README of the pipeline that can be replaced with a DOI badge. This is something we need to set up for nf-core in general, so at the moment, you need to ask someone from the core team to help get you one from [Zenodo](https://zenodo.org/). You can just add that to the README and `commit` this (and only this!) to the `master` branch directly, or ask someone to do that for you. Remember that this is optional, but it’s still considered best practice.

[11:40](https://youtu.be/1OwkTd-P5pQ?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=700) After first release, there’s a step that’s required to enable future development. So it’s a good idea after you’ve made your release on the `master` branch to do a bump version on your `dev` branch, for example `1.0.1dev`. Remember to use the suffix `dev`. You can also update the CHANGELOG file with a new section for `1.0.1dev`. If you are working on your code in your own fork in the main repository of nf-core, you need to open a `pull request` to `dev` in your pipeline from your fork with all the changes I just mentioned.

[12:27](https://youtu.be/1OwkTd-P5pQ?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=747) There are a couple of extra things that you can follow up on. You likely want to do subsequent releases of your pipeline, and there’s a [tips and tricks section](https://nf-co.re/docs/contributing/adding_pipelines#subsequent-releases) on the nf-core website, that’s specifically targeted at this. There’s also an [extra guide](https://nf-co.re/docs/contributing/adding_pipelines#adding-new-pipeline-features-to-existing-pipelines) on adding features to existing pipelines.

[13:10](https://youtu.be/1OwkTd-P5pQ?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=790) So most of what I’ve talked about is based on two how-to guidelines on the nf-core homepage (<https://nf-co.re/docs/contributing/adding_pipelines#making-the-first-release> and <https://nf-co.re/docs/contributing/release_checklist>), so do check those out.

</details>
