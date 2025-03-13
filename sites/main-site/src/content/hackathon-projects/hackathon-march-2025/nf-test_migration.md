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

- [Batch 1](https://github.com/nf-core/modules/issues/7529)
- [Batch 2](https://github.com/nf-core/modules/issues/7654)

2. Move all nf-core/subworkflows to nf-test

- [Batch 1](https://github.com/nf-core/modules/issues/7575)

## More informations

- Pytest to nf-test migration [documentation](https://nf-co.re/docs/tutorials/tests_and_test_data/migrate_to_nf-test)
- [nf-test bitesize talk](https://www.youtube.com/watch?v=K9B7JRkMpQ4)
- [nf-test guidelines](https://nf-co.re/docs/tutorials/tests_and_test_data/nf-test_writing_tests)

## How to proceed

As there are a lot of modules (i.e. 162) and subworkflows (i.e. 23) to migrate we need to stay organised to avoid duplicated work. To do so, we created three issues to track every contribution and a simple procedure as follows:

- Check if no one is already assigned to the modules/subworkflow you want to work on.
- Assign yourself to the sub-issue
- Perform the migration
- When your pull request is ready, ask for a review and link it to the sub-issue
- When merged update the hackathon project
- Rinse and repeat!
