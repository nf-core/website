---
title: Enhancing nf-core/sammyseq towards the first release
category: pipelines
slack: "https://nfcore.slack.com/archives/C025DGRJCDB"
image: "/assets/images/events/2025/hackathon-march/just-waiting-for-release.jpg"
image_alt: "A skeleton sitting on a bench. The caption above: 'Just waiting for that release'. The caption below: 'anytime soon'"
leaders:
  daisymut:
    name: Margherita Mutarelli
    slack: "https://nfcore.slack.com/team/U018GB2CSGK"
  ugoiannacchero:
    name: Ugo Maria Iannacchero
    slack: "https://nfcore.slack.com/team/U073H3NAY57"
---

## Project aim

We will work on [nf-core/sammyseq](https://github.com/nf-core/sammyseq), designed for the analysis of a brand-new NGS application to analyze the chromatin state.
We want complete it by adding a new functionality: a Hi-C-like analysis of compartments, expand and improve the alignment / preprocessing part and update documentation accordingly. Any suggestions about documentation, code readability, usability and user-friendliness is warmly welcome.
Open to contributors of all experience levels, tasks for beginners will be labeled "good first issue".

## Goals

- Add compartment analysis subworkflow and components
- Add bwa-mem alignment
- Update the schema figure (metro map) to accurately represent the current analysis workflow
- Solve as many open [issues](https://github.com/nf-core/sammyseq/issues) as we can
- Check what is needed for the first release
