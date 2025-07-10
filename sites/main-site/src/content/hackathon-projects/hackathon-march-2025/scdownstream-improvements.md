---
title: Improving the scdownstream pipeline by increasing diversity of selectable tools
category: pipelines
slack: https://nfcore.slack.com/archives/C077PUX8DBP
intro_video: ""
image: "/assets/images/events/2025/hackathon-march/scdownstream_improvements.jpg"
image_alt: "Biologist running out of ideas to fix a leaking water tank smacks duct tape(scRNAseq) over the leak"
leaders:
  KurayiChawatama:
    name: Kurayi Chawatama
    slack: "https://nfcore.slack.com/team/U08CG9KE295"
  bogrum:
    name: Emre Taha Çevik
    slack: "https://nfcore.slack.com/team/U081Q2VU4GG"
---

This project aims to enhance the [nf-core/scdownstream](https://nf-co.re/scdownstream/dev/) pipeline for downstream scRNA-seq analysis.
By expanding the range of available tools, we will improve its functionality, fostering broader adoption and making the pipeline more up-to-date and appealing to the scientific community.

## Goals

- **Expanding the range of doublet detection tools to include better-performing tools, particularly scDblFinder.**
  - _For people who have experience with scRNA-seq QC in Python or R_
  - _Based on benchmark comparisons, [scDblFinder](https://github.com/plger/scDblFinder) generally outperforms the currently listed tools._

- **Expanding the range of cell type annotation tools to include other species currently not accounted for, particularly mice.**
  - _For people who have experience with mouse or other model organism's scRNA-seq data analysis_
  - _Datasets like Tabula Muris (Senis) will be explored for support integration._

- **Improving documentation and user support to promote broader adoption.**
  - _For first-time contributors interested in making the pipeline more accessible!_

_Contributors of all experience levels are welcome!_
