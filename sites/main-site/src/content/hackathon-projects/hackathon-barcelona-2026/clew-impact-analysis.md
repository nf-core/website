---
title: Clew — impact analysis over the lineage Nextflow already writes
category: tooling
location: Barcelona
slack: https://nfcore.slack.com/team/U0B8XM3MP5X
image: "/assets/images/events/2026/hackathon-barcelona/clew-lineage-thread.png"
image_alt: "Two road signs offer the two ways to handle a bad upstream input: re-run everything blindly, or compute exactly what was affected. A red thread — the clew — leads back out."
leaders:
  QuietFlare:
    name: Seema Jagadeesh
    slack: https://nfcore.slack.com/team/U0B8XM3MP5X
---

## Goal

Run Clew on real Nextflow production or research runs from participants and see how it behaves. I am also planning to add more report adapters for nf-core pipelines. I want to collect feedback on whether this model matches how labs actually work, and if it is useful for the intended users.

## Description

There are cases where an input turns out to be defective long after the run finished: a container had a bug, a sample failed QC after its data was already used somewhere downstream, or a reference genome got updated. The usual action is to rerun the current pipeline with the `-resume` command, but it doesn't tell you which older runs used the defective input, which tasks consumed its outputs, or what to delete, re-run and disclose. 

Clew tries to answer these questions with evidence: it builds a DAG from the existing data lineage store, nf-prov RO-Crates, or `work/` symlinks, and provides an impact analysis with a remediation plan. Verdicts reference the policy rule they came from and can be recomputed offline. It also follows links between pipelines, e.g. a withdrawn sample in an rnaseq run ends up in the differentialabundance run that read its count matrix.

Try it in two minutes:

```bash
pip install clew-lineage
clew demo
```

Stdlib only, AGPL: https://github.com/QuietFlare/clew

## What participants will do

The project has something for every level:

- Run the extractors on your own runs and report what breaks (no coding).
- Tell us how samples really become unrecoverable in your lab (no coding, most valuable).
- Add an adapter for a pipeline you use (40 to 60 lines of Python).
- Extend the trigger selectors.
- Take on the hard one and give the event log an identity.
