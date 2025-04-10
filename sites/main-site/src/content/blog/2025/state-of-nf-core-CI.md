---
title: "State of the nf-core CI"
subtitle: "Run, CI, Run!"
headerImage: "https://images.unsplash.com/photo-1643279614844-afb7200c0095"
headerImageAlt: "A small child running between trees in a thick forest."
pubDate: 2025-04-09T12:00:00+01:00
authors:
  - "mashehu"
label:
  - "infrastructure"
---

## Introduction

In order to scale the nf-core organisation, we make heavy use of automation.
Every repository has a `.github/workflows/` directory with a set of YAML files that configure automated
workflows using [GitHub Actions](https://github.com/features/actions).
These do all kinds of things, the most visible being to run continuous integration (CI) tests every time
you open a pull-request or push to `master`.

Continuous integration tests are critical for validation of new code and automated quality checks that
help to prevent code regressions.
These help nf-core pipelines and modules to stay standardised and high-quality, whilst avoiding
a lot of manual code checks and providing fast feedback to developers.

## The waiting game

Of course, when something works well, we tend to use it a lot.
This is certainly true with GitHub actions and nf-core.

An average pull-request in nf-core/modules can kick off up to 20 automated jobs.
PRs for nf-core pipelines can start up even more jobs.
The nf-core GitHub organisation has an allocation of 60 action runners
(this already includes a double allocation for open source projects).

This mismatch between number of jobs launched and total concurrent capacity can
very quickly create a long queue of jobs waiting to be picked up by a runner.
This is especially noticeable during peaks of community activity, such as hackathons.

One of the most popular memes during nf-core hackathons:
![Pablo Escobar from the series Narcos waiting in different positions with the caption "Waiting](../../../assets/images/blog/state-of-nf-core-CI/ci-waiting-meme.png)

Over the years we have been working on improving the situation applying different strategies:

## 1. Self-hosted runners

Thanks to a generous donation of AWS credits, we are in theory able to run these tests on self-hosted runners on AWS.

### Manual AWS instances

We first [launched instances manually, based on demand](https://github.com/nf-core/tools/issues/1940).
This was a good first step, removing a bit of the pressure during hackathons.
But this didn't scale well, and we soon reached the runner limit also outside of hackathons.

### AWS instances managed by terraform

This lead us to search for automatically scaling solutions. After a bit of trial we had a working solution using [philips-labs/terraform-aws-github-runner/](https://github.com/philips-labs/terraform-aws-github-runner).
By using the "infrastructure as code" tool called [terraform](https://www.terraform.io/) we were able to set up AWS Lambda functions, that automatically create and destroy AWS EC2 instances with predefined configurations.
These functions were triggered by a webhook to GitHub to trigger these functions.
It did the automatic scaling for us, and we didn't need to add any extra setting files in the different nf-core repositories, like other approaches would have required.
We just needed to install the specific GitHub App in the repositories and specify `runs-on: "self-hosted"{:yml}` instead of `runs-on: "ubuntu-latest"{:yml}` in the GitHub Actions workflows.
During the Barcelona hackathon in September 2024 we even managed to [add GPU-enabled runners](https://github.com/nf-core/actions-runners/pull/10), which allowed us to test and add GPU-based tools like [parabricks](https://docs.nvidia.com/clara/parabricks/latest/index.html) to nf-core/modules and nf-core pipelines (e.g. [nf-core/sarek](https://github.com/nf-core/sarek/issues/1853)).

However, the self-hosted runners solution hasn't been perfect.
We kept running into user-permission errors, because the terraform setup came only with a rootless docker config. We also found that runners sometimes didn't pick up a job from GitHub, and debugging this was very tricky. The final straw was when the modules tests kept failing due to remnants of previous test runs and resulting permissions errors.

### Moving to runsOn

As our manual runners become more problematic, we started to search for yet another solution.
This time both [Edmund](https://github.com/edmundmiller) and me stumbled upon [runsOn](https://runs-on.com/).
This service provides an AWS CloudFormation template to create a self-hosted runner in your own AWS account. So it is basically doing what we did before, but without the terraform setup and comes with really good documentation.
We switched switched the CI tests for [nf-core/modules](https://github.com/nf-core/modules/pull/7840) to this new solution a week before the March 2025 hackathon and it worked like a charm.
It was so smooth I had to check several times that the tests are actually running on the self-hosted runners.

Switching to the new solution was as simple as changing:

```diff title="nf-core/modules/.github/workflows/nf-test.yml"
- runs-on: "self-hosted"
+ runs-on:
+      - runs-on=${{ github.run_id }}
+      - runner=4cpu-linux-x64
```

And if we wanted to add GPU-enabled runners, we just specify a different runner in the workflow.
RunsOn furthermore allowed us to quickly support the initiative to add `linux/arm64` containers to nf-core/modules,
which was a great effort by other nf-core members, by adding a ARM-based runners to the tests, e.g. to the [nf-core/rnaseq pipeline tests](https://github.com/nf-core/rnaseq/pull/1530).

An additional benefit compared to the previous terraform setup is that we now use spot instances, which are both cheaper and faster to start up.

There were still some issues where the runners behaved differently compared to the GitHub-hosted runners:

- We had problems accessing some public AWS resource. We fixed this by setting:

```groovy title="nextflow.config"
aws.client.anonymous = true
```

- Dependencies were missing for `setup-apptainer` (which was fixed by the [runsOn developer](https://github.com/runs-on/runs-on/releases/tag/v2.7.0)).

But overall this new runner setup works really well.

## 2. nf-test sharding

A further approach on speeding up the CI tests was introduced in [nf-test 0.9.0](https://github.com/askimed/nf-test/releases/tag/v0.9.0).
This release added the option to split up the tests into multiple jobs using the `--shard` flag.
Instead of running all tests one after the other, nf-test will now distribute the tests into multiple jobs and run them in parallel based on given maximum number of shards.
This is especially useful for modules/subworkflows/pipelines that trigger a large number of tests.

Let's look at the tests for the [FASTQC module](https://github.com/nf-core/modules/blob/master/modules/nf-core/fastqc/tests/main.nf.test) as an example.
This module contains 12 tests, e.g. `sarscov2 single-end [fastq]`,`sarscov2 paired-end [fastq]` and their stub versions (and it is also included in 2 subworkflows, which would get tested during an update).
Without sharding, we would run all 12 tests sequentially, which can take quite a while.
With sharding, we can now run, for example, 4 shards/jobs in parallel, each with 3 tests and the tests from the subworkflows distributed over them.
This will not lower the number of needed GitHub Actions runners, but we tackled that problem with self-hosted runners.
Sharding will instead reduce the time it takes to run the tests for a single PR.

To handle this scaling, we added [a sharding step](https://github.com/nf-core/modules/blob/master/.github/actions/get-shards/action.yml) to the CI workflows, that first gets the number of triggered tests by running nf-test in dry-run-mode.
We then use this number to set the number of shards needed based on the number of tests and the `max_shard_size` parameter, which gives us a bit more control over the number of runners used and to avoid idle runners.

The disadvantage of this approach is that it is not immediately clear which tests failed during a run, one needs to instead go through the logs of the different jobs.
We tried to mitigate this by adding a summary of the tests in the GitHub Actions, but this needs some more polishing.

Thanks to the work by [@GallVp](https://github.com/GallVp), [@sateeshperi](https://github.com/sateeshperi),
[@edmundmiller](https://github.com/edmundmiller), [@adamrtalbot](https://github.com/adamrtalbot), [@mirpedrol](https://github.com/mirpedrol),
[@maxulysse](https://github.com/maxulysse), and [@mashehu](https://github.com/mashehu),
we now have this dynamic sharding step in the CI workflows for
[nf-core/modules](https://github.com/nf-core/modules/blob/ab281000d296e2a6ab4efceb14c1151bd3a326da/.github/actions/get-shards/action.yml) and in the
upcoming [pipeline template](https://github.com/nf-core/tools/blob/ae7760bcff980809f6dabdcaa96209b60a3d2d5a/nf_core/pipeline-template/.github/actions/get-shards/action.yml) for nf-core/tools 3.3.0.

### A different waiting game

Both of these strategies already removed quite a bit of waiting time and hopefully in the future we can return to another hackathon evergreen meme:

![Bernie Sanders meme with the caption: "I am once again asking for PR reviews"](../../../assets/images/blog/state-of-nf-core-CI/pr-review-bernie-meme.png)
