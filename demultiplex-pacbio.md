---
title: Adding PacBio Demultiplex Support to nf-core Pipelines
category: pipelines
slack: https://nfcore.slack.com/archives/CL88J906S
intro_video: ""
image: "/assets/images/events/2025/hackathon-march/pacbio-demultiplex.jpg"
image_alt: "Illustration of PacBio sequencing with long-read (lima) and short-read (obc2fastq) support."
leaders:
  glicht:
    name: Gabriel Lichtenstein
    slack: "https://nfcore.slack.com/team/U043QKAUR2T"
---

This project proposes the integration of PacBio demultiplex support into nf-core pipelines—addressing both long-read and short-read sequencing. We plan to use [lima](https://github.com/PacificBiosciences/barcoding) for long-reads and [obc2fastq](https://www.pacb.com/onso/software-downloads) for short-reads. Designed as an innovative hackathon project, this initiative expands nf-core's capabilities to meet the evolving needs of sequencing analyses.

## Goals

- **Integrate PacBio demultiplexing for long-reads using lima.**
  - _For users working with single-molecule sequencing data producing long-reads._
- **Integrate PacBio demultiplexing for short-reads using obc2fastq.**
  - _For users interested in processing of short-read data from the ONSO platform._
- **Ensure compatibility with nf-core and Nextflow standards.**
  - _For contributors experienced in workflow development aiming at high reliability._
- **Benchmark and validate the performance of the new functionalities during the hackathon.**
  - _For users keen on evaluating cutting-edge sequencing data processing tools._
- **Enhance documentation and user support for the new features.**
  - _For both first-time contributors and seasoned developers aiming to improve accessibility._

_We welcome contributors of all experience levels to collaborate and drive this exciting project forward during the nf-core hackathon!_
