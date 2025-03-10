---
title: Improving the References pipeline
category: pipelines
slack: "https://nfcore.slack.com/channels/references"
leaders:
  maxulysse:
    name: Maxime U Garcia
    slack: "https://nfcore.slack.com/team/UE6D8290F"
---

This project aims to enhance the:

- [nf-core/references](https://nf-co.re/references/dev) pipeline with a focus on the [1.0.0](https://github.com/nf-core/references/milestone/1) milestone.

Follow @maxulysse in a deep dive down the rabbit hole.

## Goals

- Make the schema coherent (file names)
- Complete rnaseq indices generation
- Complete sarek indices generation
- Enhance schema to handle multiple files
- Figure out how to deal with a variable base_path
- Auto-generate params.yml for backwards compatibility for old rnaseq versions
- Improve utils suborkflow usage
