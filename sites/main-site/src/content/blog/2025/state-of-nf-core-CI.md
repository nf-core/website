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

This lead us to search for automatically scaling solutions. After a bit of trial we had a working solution using [philips-labs/terraform-aws-github-runner/](https://github.com/philips-labs/terraform-aws-github-runner).
It did the automatic scaling for us, and we didn't need to add any extra setting files, just needed to install the specific GitHub App in the needed repositories and specify the self-hosted runner in the `runs-on` field in the GitHub Actions workflow.
During the Barcelona hackathon in September 2024 we even managed to add GPU-enabled runners.
But we kept running into permission errors, because the terraform setup came only with a rootless docker config. Debugging why the runners sometimes didn't pick up a job was also tricky. The final straw was when the modules tests kept failing due to remnants of previous test runs and resulting permissions errors.

So we again searched for another solution. This time both [Edmund](https://github.com/edmundmiller) and me stumbled upon [runsOn](https://runs-on.com/).
It is basically doing what we did before, without the terraform setup and comes with really, really good documentation.
We switched the [nf-core/modules](https://github.com/nf-core/modules/pull/7840) to this new solution a week before the March 2025 hackathon and it worked like a charm.
It was so smooth I had to check several times that the tests are actually running on the self-hosted runners. Switching to the new solution was really just changing:

```diff title="nf-core/modules/.github/workflows/nf-test.yml"
- runs-on: "self-hosted"
+ runs-on:
+      - runs-on=${{ github.run_id }}
+      - runner=4cpu-linux-x64
```

in the GitHub Actions workflow.
Using GPU runners was a similar change.
Because we were already at it and because of the great effort on adding `linux/arm64` containers to nf-core/modules, we also could also quickly add ARM-based CI runners to the nf-core/rnaseq pipeline tests.
An additional benefit compared to the previous terraform setup is that we now use spot instances, which are both cheaper and faster to start up.
There were still some issues where the runners behaved differently compared to the GitHub-hosted runners. We saw for example problems accessing some public AWS resources (setting `aws.client.anonymous = true{:groovy}` in the `nextflow.config` fixed this) and dependencies were missing for `setup-apptainer` (which was fixed by the [runsOn developer](https://github.com/runs-on/runs-on/releases/tag/v2.7.0)).
But overall this new runner setup works really well.

## 2. nf-test sharding

[nf-test 0.9.0](https://github.com/askimed/nf-test/releases/tag/v0.9.0) introduced the option to split up the tests into multiple jobs using the `--shard` flag.
This is especially useful for modules/subworkflows/pipelines that trigger a large number of tests.
Instead of running all tests one after the other, nf-test will now split up the tests into multiple jobs and run them in parallel based on given maximum number of shards.
We also added [a logic](https://github.com/nf-core/modules/blob/master/.github/actions/get-shards/action.yml) to automatically detect the number of shards needed based on the number of tests and the `max_shard_size` parameter.

Both of these strategies already removed quite a bit of waiting time and hopefully in the future we can return to another hackathon tradition instead:

![Waiting for PR review](../../../assets/images/blog/state-of-nf-core-CI/pr-review-waiting-meme.png)
