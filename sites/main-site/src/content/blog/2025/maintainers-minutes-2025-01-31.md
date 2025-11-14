---
title: "Maintainers Minutes: January 2025"
subtitle: "Keeping you informed of the latest maintainers discussions"
pubDate: 2025-01-31T12:00:00+01:00
headerImage: "/assets/images/blog/maintainers-minutes-2024-07-26/maintainers-wide.png"
headerImageAlt: "Cartoon yellow rubber duck with nf-core logo badge on its body with the nf-core logo."
embedHeaderImage: true
authors:
  - "maxulysse"
label:
  - "maintainers"
---

import Profile from "@components/GitHubProfilePictureExtended.astro";
import { Image } from "astro:assets";

The 'Maintainers Minutes' aims to give further insight into the workings of the [nf-core maintainers team](/governance#maintainers)
by providing brief summaries of the monthly team meetings.

## Overview

In the January meeting we discussed (amongst others) the following topics:

- [Spring cleaning](#spring-cleaning)
- [References](#references)
- [stubs](#stubs)
- [CI with latest-everything](#ci-with-latest-everything)
- [Helpdesk Americas](#helpdesk-americas)
- [Default prefixes in modules](#default-prefixes-in-modules)
- [Best practices for spring cleaning](#best-practices-for-spring-cleaning)

## Spring cleaning

We decided that the spring cleaning that we planned to do before the hackathon should be done two weeks before the hackathon.
So it will happen from March 10th to March 14th.
We decided to extend the duration of the spring cleaning to a full week to give everyone enough time to participate.
The idea is to do some actual spring cleaning in the nf-core GitHub organization.
We will go trough the common repos, update old Issues (fixed, still relevant or won't do...), close old PRs (with warning), delete branches that have already been merged...
We will check if teams (GitHub/Slack) and make sure they are up to date.
We will also look at pipelines and see if they are still active or need to be archive

## References

Maxime Garcia wanted to know about the usage of `conf/igenomes.config`, and how it was actually used in the pipelines.
Since no one is really actually using it, he will move forward with his ideas, and will continue to work on his current POC.

> Whenever Maxime says “I can do whatever I want” I get really scared.

## Stubs

We discussed stubs and how we use it.
One concern is that the Nextflow language server was complaining about unused args in the stub section of any modules
So people started to remove them, but then we have a stub block that is widely different from the script block, and no really way of knowing if the args are passed correctly within stubs.
The discussion about stubs was going a bit in circles and off topic, but we still managed to reach an agreement, and for now:

- We will echo args for now instead of removing them
- Talk with Ben, in particular in the static type PR and see if we can be smarter about it

## CI with latest-everything

We widely agreed that since the `latest-everything` Nextflow version we use in nf-core CI can sometimes break a pipeline due a change of syntax, behavior, infrastructure usage...
It is important to test with it (I do think that was one of the primary reason why we started doing it), because we want our code to be up to date with the latest update.
But we also agreed that we don't want failure in this tests to block PRs.
The end idea is that we keep the `latest-everything` tests in the CI workflow, but we don't require them to merge the PR.

## Helpdesk Americas

We discussed the helpdesk in the Americas timezone.
Since the start of the year, Florian and Lorena have mainly been on their own with very little audience.
We thought about ideas to improve it, and came up with the following:

- Consider moving the helpdesk to a later time slot so that it is outside of the busiest hours.
- Pump up on the communication, and make it more visible
  - Explaining clearly what the helpdesk is, and that there are clearly no stupid questions
  - See with the outreach team if we can get more visibility

Since we have all our sessions planed ahead, we can maybe communicate about it in a more efficient manner, as we currently rely on the people hosting the helpdesk to communicate about it.

## Default prefixes in modules

Our amazing modules aficionado Simon kindly reminded us about [this PR](https://github.com/nf-core/website/pull/2608).
It was decided that four maintainers should review it, so that we can hopefully close the discussion and advance on this topic and merge the PR.

## Best practices for spring cleaning

We discussed Spring Cleaning once more, and our relentless modules dark knight Simon shared a nice idea:

> We should check our own open PRs and close them if they are no longer relevant.

It would be nice indeed to keep track of all our open PRs, and issues and close them if we can.
That would simplify the spring cleaning, and all of our Justice League of maintainers' own lives.

## The end

As always, if you want to get involved and give your input, join the discussion on relevant PRs and slack threads!

\- :heart: from your #maintainers team!
