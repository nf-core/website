---
title: nf-core/phaseimpute updates
category: pipelines
slack: "https://nfcore.slack.com/archives/C05G7UE94CT"
intro_video: ""
image: ""
image_alt: ""
location: "Brest (GGB) and Online"
leaders:
  louislenezet:
    name: Louis Le Nézet
    slack: "https://nfcore.slack.com/team/U03UCQH8FN2"
---

This project aims to update the [nf-core/phaseimpute](https://nf-co.re/phaseimpute/dev) pipeline, by adding features and fixing bugs described in the [issues](https://github.com/nf-core/phaseimpute/issues) page.

## Goal

Update the logic and add new features to the pipeline.

## Tasks

See below a couple of potential new features to add during the hackathon.
Depending on the participation and level of the attendees, more issues could be made available for the hackathon.

### Update to topic usage and strict syntax

With the new nextflow linting tool, we aim for 0 errors and 0 warnings.

Level: _For anyone interested. Beginners welcome!_

### Improve documentation

The documentation would need some proof checking !

Level: _For anyone interested. Beginners welcome!_

### Use nf-core/subworkflow

Multiple local subworkflow should be changed to nf-core one.
Most of them are already on the nf-core/modules repository.

Level: _For people with some experience in Nextflow/nf-core pipelines._

### Add `QUILT2` imputation

Issue [#116](https://github.com/nf-core/phaseimpute/issues/116)

The `QUILT` software is now available in a new version that would need to be added.

### Add sexual chromosome imputation

Issue [#31](https://github.com/nf-core/phaseimpute/issues/31)

Update the diploidy information depending on the chromosomal position.
Set a parameter to select the start and end of the pseudo-autosomal region.

Level: _For people with experience in genomic imputation._

### Fix bugs and other issue

There is other [issues on github](https://github.com/nf-core/phaseimpute/issues), feel free to check them out or add yours if you encounter any.
