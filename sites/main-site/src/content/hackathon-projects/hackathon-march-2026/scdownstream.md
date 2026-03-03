---
title: Improving the scdownstream pipeline by increasing diversity of selectable tools
category: pipelines
slack: https://nfcore.slack.com/archives/C077PUX8DBP
intro_video: ""
image: "/assets/images/events/2026/hackathon-march/2-3-years-dagestan-meme.jpg"
image_alt: "Biologist running out of ideas to fix a leaking water tank smacks duct tape(scRNAseq) over the leak"
leaders:
  KurayiChawatama:
    name: Kurayi Chawatama
    slack: "https://nfcore.slack.com/team/U08CG9KE295"
  bogrum:
    name: Emre Taha Çevik
    slack: "https://nfcore.slack.com/team/U081Q2VU4GG"
---

## Goals

This project aims to further enhance the [nf-core/scdownstream](https://nf-co.re/scdownstream/dev/) pipeline by implementing new analytical features and modernizing existing components.

By addressing open feature requests and technical debt, we aim to improve reproducibility, usability, and alignment with current single-cell best practices.

_Contributors of intermediate and advanced experience levels are especially encouraged to participate._

## Tasks

- **Implement cell cycle analysis** [(issue #170)](https://github.com/nf-core/scdownstream/issues/170)
  - Integrate standardized cell cycle scoring into the pipeline
  - Ensure compatibility with existing downstream modules

- **Implement `scDblFinder` as a module** [(issue #144)](https://github.com/nf-core/scdownstream/issues/144)
  - Add `scDblFinder` as an nf-core compliant module
  - Ensure proper container/Conda environment support

- **Add sample sex prediction functionality** [(issue #171)](https://github.com/nf-core/scdownstream/issues/171)
  - Implement automated sex prediction based on expression markers
  - Integrate outputs into summary reports

- **Migrate legacy steps to use `AnndataR`** [(issue #207)](https://github.com/nf-core/scdownstream/issues/207)
  - Replace older data handling steps with `AnndataR`-based workflows
  - Improve maintainability and future compatibility of the pipeline

- **Update MAD-based dynamic filtering option PR to comply with latest nf-core template** [(PR #204)](https://github.com/nf-core/scdownstream/pull/204)
  - Refactor open pull request to align with current nf-core pipeline standards
  - Ensure template compliance and successful CI checks
