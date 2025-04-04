---
title: "State of the nf-core CI"
subtitle: "Run, CI, Run!"
headerImage: "https://images.unsplash.com/photo-1643279614844-afb7200c0095"
headerImageAlt: "A small child running between trees in a thick forest."
pubDate: 2025-04-04T12:00:00+01:00
authors:
  - "mashehu"
label:
  - "infrastructure"
---

## The waiting game

One of the most popular memes during nf-core hackathons:
![Pablo Escobar from the series Narcos waiting in different positions with the caption "Waiting](../../../assets/images/blog/state-of-nf-core-CI/ci-waiting-meme.png)

This is due to the fact that the nf-core GitHub Organisation has an allocation of 60 runners (which already includes a double allocation for open source projects).

An average PR in nf-core/modules can kick up to 20 jobs off, while pipelines like nf-core/viralrecon start 178 jobs per PR.
This very quickly creates a long queue of jobs waiting to be picked up by a runner.

Over the years we have been working on improving the situation applying different strategies:

## 1. Self-hosted runners

Thanks to a generous donation of AWS credits, we are in theory able to run these tests on self-hosted runners on AWS.

We first [launched instances manually, based on demand](https://github.com/nf-core/tools/issues/1940).
This was a good first step, removing a bit of the pressure during hackathons.
But this didn't scale well, and we soon reached the runner limit also outside of hackathons.

This lead us to search for automatically scaling solutions. After a bit of trial we had a working solution using [philips-labs/terraform-aws-github-runner/](https://github.com/philips-labs/terraform-aws-github-runner). By using the "infrastructure as code" tool called [terraform](https://www.terraform.io/) we were able to set up AWS Lambda functions, that automatically create and destroy AWS EC2 instances with predefined configurations.
These functions were triggered by a webhook to GitHub to trigger these functions.
It did the automatic scaling for us, and we didn't need to add any extra setting files in the different nf-core repositories, like other approaches would have required. We just needed to install the specific GitHub App in the repositories and specify `runs-on: "self-hosted"` instead of `runs-on: "ubuntu-latest"` in the GitHub Actions workflows.

During the Barcelona hackathon in September 2024 we even managed to [add GPU-enabled runners](https://github.com/nf-core/actions-runners/pull/10), which allowed us to test and add GPU-based tools like [parabricks](https://docs.nvidia.com/clara/parabricks/latest/index.html) to nf-core/modules and nf-core pipelines (e.g. [nf-core/sarek](https://github.com/nf-core/sarek/issues/1853)).
But we kept running into permission errors, because the terraform setup came only with a rootless docker config. Debugging why the runners sometimes didn't pick up a job was also tricky. The final straw was when the modules tests kept failing due to remnants of previous test runs and resulting permissions errors.

So we started to search for yet another solution. This time both [Edmund](https://github.com/edmundmiller) and me stumbled upon [runsOn](https://runs-on.com/).
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
Because we were already at it and because of the great effort by other nf-core members on adding `linux/arm64` containers to nf-core/modules, we also could quickly add ARM-based CI runners to the [nf-core/rnaseq pipeline tests](https://github.com/nf-core/rnaseq/pull/1530).

An additional benefit compared to the previous terraform setup is that we now use spot instances, which are both cheaper and faster to start up.

There were still some issues where the runners behaved differently compared to the GitHub-hosted runners.

- We had problems accessing some public AWS resources (setting `aws.client.anonymous = true{:groovy}` in the `nextflow.config` fixed this)
- Dependencies were missing for `setup-apptainer` (which was fixed by the [runsOn developer](https://github.com/runs-on/runs-on/releases/tag/v2.7.0)).

But overall this new runner setup works really well.

## 2. nf-test sharding

A further approach on speeding up the CI tests was introduced in [nf-test 0.9.0](https://github.com/askimed/nf-test/releases/tag/v0.9.0).
This release added the option to split up the tests into multiple jobs using the `--shard` flag.
Instead of running all tests one after the other, nf-test will now distribute the tests into multiple jobs and run them in parallel based on given maximum number of shards.
This is especially useful for modules/subworkflows/pipelines that trigger a large number of tests.

Let's look at the tests for the [FASTQC module](https://github.com/nf-core/modules/blob/master/modules/nf-core/fastqc/tests/main.nf.test) as an example. This module contains 12 tests, e.g. `sarscov2 single-end [fastq]`,`sarscov2 paired-end [fastq]` and their stub versions (and it is also included in 2 subworkflows, which would get tested during an update).
With sharding enabled, we can now run, for example, 4 shards/jobs in parallel, each with 3 tests and the tests from the subworkflows distributed over them.

To handle this scaling, we added [a sharding step](https://github.com/nf-core/modules/blob/master/.github/actions/get-shards/action.yml) to the CI workflows, that first gets the number of triggered tests by running nf-test in dry-run-mode.
We then use this number to set the number of shards needed based on the number of tests and the `max_shard_size` parameter, which gives us a bit more control over the number of runners used and to avoid idle runners.

The disadvantage of this approach is that it is not immediately clear which tests failed during a run, one needs to instead go through the logs of the different jobs.
We tried to mitigate this by adding a summary of the tests in the GitHub Actions, but this needs some more polishing.

### A different waiting game

Both of these strategies already removed quite a bit of waiting time and hopefully in the future we can return to another hackathon tradition instead:

![Bernie Sanders meme with the caption: "I am once again asking for PR reviews"](../../../assets/images/blog/state-of-nf-core-CI/pr-review-bernie-meme.png)
