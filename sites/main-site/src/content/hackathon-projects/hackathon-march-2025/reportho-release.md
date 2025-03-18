---
title: Finalizing nf-core/reportho 1.1.0
category: pipelines
slack: "https://nfcore.slack.com/archives/C06R4FC4C1E"
intro_video: ""
image: "/assets/images/events/2025/hackathon-march/reportho_world.jpeg"
image_alt: "A drawing of a beautiful futuristic city, captioned with the words 'the world if reportho'"
leaders:
  itrujnara:
    name: Igor Trujnara
    slack: "https://nfcore.slack.com/archives/D03DEGFQ0AG"
---

This project aims to close open PRs and issues of the [nf-core/reportho](https://nf-co.re/reportho/dev) pipeline to finalize release 1.1.0.

## Goals

- **Update custom report**

We need to update the custom report to include the sequence-based merging step results (if run) and exclude alignment/phylogeny results.
The code for this part is partially in a [separate repo](https://github.com/itrujnara/orthologs-report/tree/main), so contact me before you work on this.

- **Update MultiQC report**

We need to make sure that the MultiQC report contains all required information and no more (i.e. no references to alignment or phylogeny).

- **Update metro map**

The workflow has changed since 1.0. We need to update the metro map to reflect that.

- **Verify documentation**

We need to verify that all documentation files are up to date. I have already made some changes, but there might be something missing.

- Finally, **submit for review**

## Note

I will not be available on Monday (March 24). Feel free to work on your forks and submit PRs during that day. I will review all changes on Tuesday morning (CET).

_We welcome contributors of all experience levels. Feel free to get in touch on the [Slack channel](https://nfcore.slack.com/archives/C06R4FC4C1E) if you have any questions or comments._
