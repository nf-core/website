---
title: nf-core/sarek language updates
category: pipelines
slack: https://nfcore.slack.com/archives/C05V9FRJYMV
intro_video: ""
location: online
image: "/assets/images/events/2026/hackathon-march/nfcore-sarek-icon.png"
image_alt: "nf-core/sarek logo"
leaders:
  FriederikeHanssen:
    name: Friederike Hanssen
    slack: https://nfcore.slack.com/team/UPC8CHDKQ
---

This project aims to update [nf-core/sarek](https://nf-co.re/sarek/dev) with newer Nextflow features.

## Goal

Adopt topic channels and finalize the strict syntax upgrade.

## Tasks

### Implement topic channels for versions and MultiQC

Replace explicit channel wiring for version collection and MultiQC inputs with Nextflow topic channels.

### Finalize strict syntax upgrade

Make sure the pipeline runs cleanly with `nextflow.enable.strict = true` and resolve any remaining warnings.

:::note
Keep an eye on this page for more updates.
:::
