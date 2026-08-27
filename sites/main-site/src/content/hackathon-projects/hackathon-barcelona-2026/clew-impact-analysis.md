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

Test Clew's impact analysis against participants' real Nextflow runs, add adapters for more nf-core pipelines, and collect practitioner corrections to the remediation model, so the tool leaves the hackathon proven against reality rather than fixtures.

## Description

When something upstream turns out to be bad (a buggy container, a failed QC, a stale reference, a withdrawn sample), `-resume` fixes the run in front of you, but nothing answers which past runs used the bad input, which other runs consumed their outputs, and what has to be deleted, re-run, or disclosed now.

Clew computes that. It reads the lineage Nextflow already writes (the native lineage store on 25.04+, nf-prov RO-Crates, or plain `work/` symlinks on any version), merges runs into one graph, and produces a remediation plan where every verdict cites a rule and replays offline, so the output is something you can hand to an auditor. It traces across pipelines too: a withdrawn sample in an rnaseq run reaches the differentialabundance run that consumed its count matrix.

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
