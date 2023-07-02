---
title: 'Bytesize: Workflow safety and immutable objects'
subtitle: Rob Syme, Seqera Labs
type: talk
start_date: '2023-05-23'
start_time: '13:00 CET'
end_date: '2023-05-23'
end_time: '13:30 CET'
location_url:
  - https://seqera-io.zoom.us/j/89853414209?pwd=V2NoK1lvN0FGMDRhNW4rbXlhWVRvQT09
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: Workflow safety and immutable objects

This week, Rob Syme ([@robsyme](https://github.com/robsyme)) will talk about how to avoid introducing subtle concurrency bugs in your Nextflow workflows by safely modifying objects (or more specifically, by _not_ modifying objects) when writing Groovy closures.
