---
title: 'Bytesize: Get your containers with "Nextflow inspect"'
subtitle: Phil Ewels, Seqera
type: talk
startDate: "2024-05-28"
startTime: "13:00+02:00"
endDate: "2024-05-28"
endTime: "13:30+02:00"
locations:
  - links:
      - https://kth-se.zoom.us/j/68390542812
---

With the advent of DSL2, nf-core pipelines began to use separate containers for every process.
Now, many nf-core pipelines can have upwards of 50 different software containers in a single run.
Traditionally, finding the URIs for each of these containers was a difficult task, with several people developing tooling based around regular expressions.

In Nextflow `23.09.0-edge` a new subcommand was added called `nextflow inspect`.
This command retrieves all of the containers that will be used by a pipeline run and returns them in Nextflow config syntax.

In this talk, Phil Ewels ([@ewels](https://github.com/ewels/)) will guides us through the process of fetching container URIs using `nextflow inspect`
and show how this can be used to customise containers in a pipeline. This may be of particular use to people working in environments that cannot
access our public container registries or who may need additional steps in their container verification steps.

**Unfortunately, the recordings for this bytesize was corrupted. Some of the content, however, was repeated in the bytesize [Bytesize: Explaining Wave containers](https://nf-co.re/events/2024/bytesize_using_wave).**
