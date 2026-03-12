---
title: Improving and fixing issues in nf-core/hic
category: pipelines
location: "ZS Copenhagen"
slack: "https://nfcore.slack.com/archives/CJ3UVAZFT"
leaders:
  samuelruizperez:
    name: Samuel Ruiz-Pérez
    slack: "https://nfcore.slack.com/team/U05H41A4BDJ"
    github: "https://github.com/samuelruizperez"
---

This project aims to address open issues and improve the [nf-core/hic](https://nf-co.re/hic/) pipeline for Hi-C data analysis. We will address the milestones for the [release of pipeline version 2.2.0](https://github.com/nf-core/hic/milestone/1).

The nf-core/hic pipeline processes Hi-C data to analyze chromatin 3D organization. With 42+ open issues covering bugs, enhancements, and feature requests, this hackathon project provides opportunities to contribute at all levels.

We welcome contributors of all experience levels!

## Goals

### Bug fixes

- Fix parameter type errors ([#229](https://github.com/nf-core/hic/issues/229): `--res_zoomify` type validation)
- Resolve MultiQC reporting issues ([#221](https://github.com/nf-core/hic/issues/221): not all samples reported)
- Fix compressed genome file handling ([#215](https://github.com/nf-core/hic/issues/215): GET_RESTRICTION_FRAGMENTS)
- Address bowtie2 version compatibility ([#208](https://github.com/nf-core/hic/issues/208))
- Resolve COOLER_MAKEBINS process errors ([#218](https://github.com/nf-core/hic/issues/218))
- Fix FastQ splitting errors ([#198](https://github.com/nf-core/hic/issues/198))
- And more...

### Enhancements

- Improve Arima Hi-C protocol support ([#225](https://github.com/nf-core/hic/issues/225))
- Add loop calling tools ([#189](https://github.com/nf-core/hic/issues/189))
- Improve MultiQC report ([#190](https://github.com/nf-core/hic/issues/190))
- Compress output files (.allValidPairs, .txt) ([#220](https://github.com/nf-core/hic/issues/220))
- Implement replicate merging functionality ([#223](https://github.com/nf-core/hic/issues/223))
- Add pairtools protocol support ([#162](https://github.com/nf-core/hic/issues/162))

### Documentation and usability

- Update documentation to reflect current pipeline state
- Improve error messages
- Add examples for common use cases

## Recommended preparation

To help you get the most out of this project, we suggest completing these Nextflow training modules beforehand:

- [Hello Nextflow](https://training.nextflow.io/latest/hello_nextflow/)
- [Nextflow Run](https://training.nextflow.io/latest/nextflow_run/)
- [Hello nf-core](https://training.nextflow.io/latest/hello_nf-core/)

## Resources

- [nf-core/hic repository](https://github.com/nf-core/hic)
- [Open issues](https://github.com/nf-core/hic/issues)
- [Pipeline documentation](https://nf-co.re/hic/)
