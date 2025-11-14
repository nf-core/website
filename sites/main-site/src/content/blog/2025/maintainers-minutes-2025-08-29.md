---
title: "Maintainers Minutes: August 2025"
subtitle: "Keeping you informed of the latest maintainers discussions"
pubDate: 2025-08-31T10:00:00+01:00
headerImage: "/assets/images/blog/maintainers-minutes-2024-07-26/maintainers-wide.png"
headerImageAlt: "Cartoon yellow rubber duck with nf-core logo badge on its body with the nf-core logo."
embedHeaderImage: true
authors:
  - "maxulysse"
label:
  - "maintainers"
---

The 'Maintainers Minutes' aims to give further insight into the workings of the [nf-core maintainers team](/governance#maintainers) by providing brief summaries of the monthly team meetings.

## Overview

Summer was actually not over as we only were a small group on that Friday.

We talked briefly about the aftermath of the previous meeting regarding the test datasets.
Edmund (@edmundmiller) reported that Cloudflare should be cheap, but we would need a wider room of people to discuss the pros and cons of the two leading options of AWS S3 vs. Cloudflare.

And then we talked about plans for migrating the code base of modules.
Upcoming in the nearish™️ future, we are expecting migrations including:

1. Linting and automatic formatting of all module code via the Nextflow language server
2. Adoption of 'topic' channels, in particular for `versions` reporting
3. Move towards Seqera containers

In all cases the various new functionality is not 💯 ready for roll out yet, for example with Harshil Alignment not completely working in the Nextflow formatter, or our experimental implementations not yet being fully tested, such as with topics being implemented in the nf-core-utils plugin.

Therefore there is no new news, but watch this space!

At the end, Maxime as his usual self, showcased what he has been playing around with lately, and demoed nf-core-utils and the particular nf-test setup in sarek.

## The end

We will continue to work on these topics in the coming months.

As always, we are looking for community involvement, so we encourage anyone to join the discussion in relevant PRs and Slack channels and threads!

\- :heart: from your #maintainers team!
