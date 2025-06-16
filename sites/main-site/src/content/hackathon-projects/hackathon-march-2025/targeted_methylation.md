---
title: Adding targeted analysis to nf-core/methylseq pipeline
category: pipelines
slack: https://nfcore.slack.com/archives/CP3RJSMF0
intro_video: ""
image: ""
image_alt: ""
leaders:
  luigilamparelli:
    name: Luigi Lamparelli
    slack: https://nfcore.slack.com/team/U089MLSGA2Y
  dcarrillox:
    name: Daniel Carrillo Bautista
    slack: https://nfcore.slack.com/team/U08AYP59MD5
---

This project focuses on developing a targeted analysis for methylseq data in the [nf-ccore/methylseq](https://github.com/nf-core/methylseq) pipeline.
We’re planning to add the following features:

- `--targeted` parameter: Currently, all detected methylation signals are included in the final results. However, some enrichment protocols target specific genomic regions. To address this, we propose adding a `--targeted` parameter that allows filtering results to retain only on-target methylation signals. The input would be a BED file specifying the target regions.
- Enrichment metrics: Building on this idea, we propose adding a job to calculate on-target vs. off-target fractions of aligned reads, along with other enrichment metrics. A tool like `CollectHsMetrics` from Picard could be used for this purpose. There’s already an open [PR](https://github.com/nf-core/methylseq/pull/313) related to this, which we’d be happy to contribute to.

## Goals

- **Enrichment filtering results.**
  - _For people interested in Nextflow module development and targeted methylation sequencing analysis._
- **Ensuring smooth integration with existing pipeline components.**
  - _For contributors with experience in nf-core module structure and dependencies._
- **Create enrichment metrics using Picard `CollectHsMetrics`.**
  - _For users that would like to integrate an existing nf-core module into a new pipeline._

_We welcome contributors of all experience levels._
