---
title: nf-test migration
category: components
slack: "https://nfcore.slack.com/archives/C049MBCEW06"
intro_video: "https://www.youtube.com/watch?v=MuwluoKFnwM"
image: "/assets/images/events/2025/hackathon-march/nf-test.png"
image_alt: nf-test everything
leaders:
  louislenezet:
    name: Louis Le Nézet
    slack: "https://nfcore.slack.com/team/U03UCQH8FN2"
---

This project aims to finalise the transition of unit tests for all nf-core modules and sub-workflows to nf-test.
It's a good starting point for first timer to contribute and review Pull-Request.

## Goal

1. Move all nf-core/modules to nf-test
2. Move all nf-core/subworkflows to nf-test

## More informations

Pytest to nf-test migration: https://nf-co.re/docs/tutorials/tests_and_test_data/migrate_to_nf-test

## How to proceed

As there is a lot of modules (i.e. 162) and subworkflows (i.e. 23) to migrate we need to stay organise to avoid duplicated work. To do so, we created a [project list](https://hackmd.io/bZ1dM18_TMqLPxFxg1NjGQ) and a simple procedure as follow:

- Check if no one is already assigned to the modules/subworkflow you want to work on.
- Add your github username to the hackmd
- Perform the migration
- When your pull request is ready, ask for a review and add the PR number to the hackmd
- When merged add the date to the hackmd
- Rinse and repeat !
