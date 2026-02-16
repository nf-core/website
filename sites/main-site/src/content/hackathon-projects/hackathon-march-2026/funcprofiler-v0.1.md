---
title: Preparing funcprofiler for first release
category: community
slack: https://nfcore.slack.com/archives/C09A1GY7J4B
# intro_video: ""
location: online/NC
leaders:
  nickp60:
    name: Nick Waters
    slack: https://nfcore.slack.com/team/U08V65MQSV9
---

We aim to release v0.1 of the read-based metagenome functional profiling `funcprofiler` pipeline.

## Goal

We aim to address some of remaining TODOs needed for submitting `funcprofiler` for review with the nf-core team.

## Description

Currently funcprofiler can successfully run `fmhfunprofiler`, `HUMAnNv3`, `HUMAnNv4`. `RGI-BWT` is in the works, pending module acceptance.
There are a number of tasks to address before the pipeline can be submitted for review.
Additionally, if users have specific tools they wish to add, this would be a great time to contribute!

## Tasks

Participating in this team could include the following tasks:

- making mini, nf-test-dataset compatible databases for tools
- submitting local modules to nf-core/modules
- implementing nf-tests
- writing documentation
- style, linting, etc
