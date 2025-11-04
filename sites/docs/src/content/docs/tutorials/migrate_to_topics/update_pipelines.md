---
title: Update pipelines to use topic channels
description: How to update pipelines to use the new topic channels feature
---

Adapting pipelines to use topic channels for version outputs is a straightforward process explained in the steps below:

1. Update the template of the pipeline to the latest version that includes support for topic channels. This is supported starting from version 3.5.0.

1. Pull the latest changes made to modules using `nf-core modules update`.

1. Look for issues in the pipeline code where `versions` is an invalid output of a process. Remove these outputs since these are now handled by the topic channels. The errors will most likely look like this:

```
ERROR ~ No such variable: Exception evaluating property 'versions' for nextflow.script.ChannelOut, Reason: groovy.lang.MissingPropertyException: No such property: versions for class: groovyx.gpars.dataflow.DataflowBroadcast
```

Example: We just updated the `samtools/sort` module to use topic channels. The pipeline will most likely have something like this in one of the (sub)workflows:

```nextflow
SAMTOOLS_SORT(input)
ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
```

The above error can be fixed by removing the line that references `SAMTOOLS_SORT.out.versions`, as shown below:

```nextflow
SAMTOOLS_SORT(input)
```
