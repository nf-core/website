---
title: nf-core/proteinfamilies updates
category: pipelines
slack: "https://nfcore.slack.com/archives/C07URFZF074"
location: Wellcome Genome Campus, Cambridge, Hinxton, UK
image: ""
image_alt: ""
leaders:
  vagkaratzas:
    name: Evangelos Karatzas
    slack: "https://nfcore.slack.com/team/U05LNHCFLCW"
---

This project aims to update the [nf-core/proteinfamilies](https://nf-co.re/proteinfamilies/dev) pipeline, by adding features described in the [issues](https://github.com/nf-core/proteinfamilies/issues) page.

## Goal

Update the logic and add new features to the pipeline.

## Tasks

See below a couple of potential new features to add during the hackathon.
Depending on the participation and level of the attendees, more issues could be made available for the hackathon.

### Update the update mechanism

Issue [#122](https://github.com/nf-core/proteinfamilies/issues/122)

Update the update mechanism of the pipeline to also output updated families `.faa` files,
similarly to the basic execution mode.

Level: _For anyone interested. Beginners welcome!_

### Update and add the trimal module to the pipeline

Issue [#62](https://github.com/nf-core/proteinfamilies/issues/62)

Update the `trimal` nf-core module and then use it in the `nf-core/proteinfamilies` pipeline.
There can be two separate tasks that use this module:

1. As a `clipkit` alternative, to trim the gappy ends of multiple sequence alignments, that is also possibly faster (to benchmark).
2. As an alternative for the strict `mmseqs` clustering, for sequence similarity removal within families (to benchmark).

Level: _For people with some experience in Nextflow/nf-core pipelines._
