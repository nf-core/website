---
title: Update pipelines to use topic channels
description: How to update pipelines to use the new topic channels feature
---

Adapting pipelines to use topic channels for version outputs is a straightforward process explained in the steps below:

1.  Update the template of the pipeline to the latest version that includes support for topic channels (nf-core tools >=3.5.0).

    The main change in the template are the following lines:

    ```groovy title="$PIPELINE_NAME.nf"
      def topic_versions = channel.topic("versions")// [!code ++]
          .distinct()// [!code ++]
          .branch { entry ->// [!code ++]
              versions_file: entry instanceof Path// [!code ++]
              versions_tuple: true// [!code ++]
          }// [!code ++]
    // [!code ++]
      def topic_versions_string = topic_versions.versions_tuple// [!code ++]
          .map { process, tool, version ->// [!code ++]
              [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]// [!code ++]
          }// [!code ++]
          .groupTuple(by:0)// [!code ++]
          .map { process, tool_versions ->// [!code ++]
              tool_versions.unique().sort()// [!code ++]
              "${process}:\n${tool_versions.join('\n')}"// [!code ++]
          }// [!code ++]
    // [!code ++]
      ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))// [!code ++]
          .mix(topic_versions_string)// [!code ++]
      ch_collated_versions = softwareVersionsToYAML(ch_versions)// [!code --]
    ```

1.  Pull the latest changes made to modules using `nf-core modules update`.

1.  Look for issues in the pipeline code where `versions` is an invalid output of a process.
    Remove these outputs since these are now handled by the topic channels. The errors will most likely look like this:

    ```console title="nextflow.log"
    ERROR ~ No such variable: Exception evaluating property 'versions' for nextflow.script.ChannelOut,
    Reason: groovy.lang.MissingPropertyException:
    No such property: versions for class: groovyx.gpars.dataflow.DataflowBroadcast
    ```

    Example: We just updated the `samtools/sort` module to use topic channels. The pipeline will most likely have something like this in one of the (sub)workflows:

    ```groovy title="main.nf"
    SAMTOOLS_SORT(input)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
    ```

    The above error can be fixed by removing the line that references `SAMTOOLS_SORT.out.versions`, as shown below:

    ```groovy title="main.nf"
    SAMTOOLS_SORT(input)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions) // [!code --]
    ```
