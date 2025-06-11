---
title: Improving the Nextflow CO2 plugin
category: tooling
slack: https://nfcore.slack.com/archives/C060UBAQ5QF
image: /assets/images/events/2025/hackathon-march/globalcore-stripes.png
image_alt: The average annual global temperature over the years 1850-2017, known as the 'warming stripes' figure from the [climate lab book](https://www.climate-lab-book.ac.uk/2018/warming-stripes/) website
leaders:
  skrakau:
    name: Sabrina Krakau
    slack: https://nfcore.slack.com/team/UMLKFJ264
  JosuaCarl:
    name: Josua Carl
    slack: https://nfcore.slack.com/archives/D08JJNCQ4SJ
  nadnein:
    name: Nadja Volkmann
    slack: https://nfcore.slack.com/archives/D08KCG9SRSL
---

# nf-co2footprint plugin project

## Resources

- [Website](https://nextflow-io.github.io/nf-co2footprint/)
- [Github](https://github.com/nextflow-io/nf-co2footprint)
- [Slack](https://nextflow.slack.com/archives/C060UBAQ5QF)

## Summary

Help us to solve problems we encountered during the path to the first release of a CO2 estimation plugin for Nextflow runs.

## Goals

- Address open issues in of the current plugin development process
- Brainstorm for potential features

## Tasks

- Add information on potential CO2 footprint optimization
  - [Issue](https://github.com/nextflow-io/nf-co2footprint/issues/110)
- Accumulate CO2 Emission of several runs from cache
  - [Issue](https://github.com/nextflow-io/nf-co2footprint/issues/62)
  - One value for included cached footprint, one for real execution value
- `cpu_model` comparison to vendor info
  - [Issue](https://github.com/nextflow-io/nf-co2footprint/issues/61)
- Compare the Plugin results to the AWS Report
- Run plugin on AWS; make it applicable for AWS
  - Info about server and server location accessible?
  - Does AWS use special AWS servers?
