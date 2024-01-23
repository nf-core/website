---
title: 'Bytesize 6: All about modules'
subtitle: Friederike Hanssen / Kevin Menden - QBiC Tübingen, Germany
type: talk
startDate: '2021-03-09'
startTime: '13:00+01:00'
endDate: '2021-03-09'
endTime: '13:30+01:00'
youtube_embed:
  - https://www.youtube.com/embed/tWvou0xj9wA
  - https://www.youtube.com/embed/Wc4A2tQ6WWY
locationURL:
  - https://youtu.be/tWvou0xj9wA
  - https://youtu.be/Wc4A2tQ6WWY
  - https://www.bilibili.com/video/BV1nN411Q7Ex
  - https://www.bilibili.com/video/BV1bz4y117me
  - https://doi.org/10.6084/m9.figshare.14185736.v1
  - https://doi.org/10.6084/m9.figshare.14185745.v1
---

# nf-core/bytesize

Join us for an episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 6: DSL2

### DSL2 - Using modules in a pipeline

This week, Friederike Hanssen ([@friederikehanssen](http://github.com/friederikehanssen/)) will present: _**DSL2 - Using modules in a pipeline.**_

The first talk will describe how to use DSL2 modules when writing an nf-core pipeline - both custom (`local`) modules, centralised from nf-core/modules and subworkflow files.

<details markdown="1"><summary>Video transcription</summary>

:::note
This text has been edited to make it more suitable for reading.
:::

#### DSL2 - Using modules in a pipeline

Thank you very much and welcome to the bytesize talk on how to use modules in a pipeline.

I would just like to briefly recap what a module is. It is an atomic process that cannot be reduced any further and usually contains a single software tool like FastQC for example and it can be used within a pipeline and also shared between different pipelines. To make use of the sharing, there is an nf-core modules repository on GitHub where you can find many of these modules already.

So to make use of this modules repository, there’s a new sub-command in nf-core/tools and for a brief recap, you can install it with `pip install nf-core` and `conda install nf-core` and if you then run `nf-core modules`, you get a list of sub-commands that you can use to interact with this repository.

[0:47](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=51) With one of the sub-commands, and with --help, you get some instructions on how to use them, so over the next couple of slides, I will briefly introduce you to some of the sub-commands that could be helpful for using modules in the pipeline.

[1:30](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=90) One of the first things you might want to do is to try and find a module that you could use, and you can use the `nf-core modules list` for that.

That will just print out all the modules currently available and they are subdivided into tool and sub-tool because many tools like `samtools` or `bcftools` have these sub-commands that are then their own module basically.

There are other tools like `fastqc` that don’t have any sub-commands. They are just a tool.

[1:54](https://youtu.be/tWvou0xj9wAlist=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=118) So if you look through the list and find a module that you would like to install, you can run `nf-core modules install` and then your pipeline directory and the tool name, and that will install the entire module for you without having to do anything else.

[2:12](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=132) It then looks like this, so on the left hand side, there is a completely fresh new pipeline that I created with the template and then I ran `nf-core modules install` and in the green part you can see where this module ends up.

The module basically creates a subfolder called nf-core software where you can find the FASTQC module, and with the three files i.e. the `functions.nf`, `main.nf` and `meta.yml`. In the `meta.yml`, you can find documentation such as what type of input and output this module takes, who wrote it and these sort of things.

In the main nf, that’s where the extra magic happens, where the fastqc is run, and in the `functions.nf` there are some helper functions that are needed.

[3:03](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=183) I guess the opposing step to that is how to remove a module that you no longer want to use or that you erroneously installed.

For that you run `nf core modules remove` and that removes this entire subfolder that had fastqc with the functions main and meta file in it. You don’t have to do anything else there.

[3:25](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=206) Then sometimes the module will get updated, maybe the software was updated and this was propagated already to the module, and you would like to use this update. So how do we ensure that first of all we have the newest version, and then how to update it?

So to check this module, you can run `nf core modules lint` on the directory to check all modules, or on a specific module like fastqc.

It will check among other things whether or not you have missed any changes. Then to update it, you currently have to remove it and reinstall it.

But for the future, Kevin is working on creating an ´update´ subcommand that you can then use.

[4:13](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=253) OK, so last but not least, what do you do if the software you’re looking for in not in nf-core/modules?

You basically have two options, the first is to add it to nf-core/modules, and the second is to create a local module. So to add it to nf-core/modules, this is usually really helpful if others use this pipeline or software as well.

However, if you’re unsure about this, you can always check the issues to see if someone else has already started working on it. Alternatively, ask in the #modules channel in Slack. Kevin will also talk about how to do this in his presentation.

[4:57](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=297) The other option is to create a local module, so this is useful for software that is specific for your pipeline.

You can run `nf core modules create`, which will create a local module for you.

If you look in this box, you can see the nf-core/modules folder that we’ve already seen before, and the local subfolder in which your tools will then live. A lot of the things that Harshil and Kevin covered in their talks is relevant here, and Kevin’s talk is more tailored towards modules, but a lot of the functionality will be similar there.

[5:36](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=336) OK, so this `nf core modules create` gives you this file here with many to-do statements and little help messages, and you can start filling out and try to get your tool to run here.

Hopefully the to-do statements help you figure out what exactly you will need to do in each step.

[5:58](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=358) So now that we have either local or nf-core/modules, we want to start writing actual workflows and pipelines.

There are two different types; the sub-workflows that are chains of multiple modules with some sort of higher-level functionality like all the qc tools that will be run on fastqc, and then the actual workflows, which are end-to-end pipelines. Then there are the DSL1 that we’ve known as large monolithic scripts and DSL2 that is a combination of modules and sub-workflows, and this is the really taking one input and producing a final output.

[6:35](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=395) On the right hand side, I have visualised the file structure of the ´viralrecon´ pipeline, and for modules we have the local and nf-core ones that I showed you before. For the sub-workflows, it will be a similar structure because some of these might be relevant to many pipelines such as qc.

Then we have the workflows and in this case, I think are hardware separated by the input data types; for Illumina and Nanopore data there are different types of workflows, which will then be called from the main nf. These workflows consist of sub-workflows and modules.

[7:14](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=434) So to use a module, you have to install the module or create a local one and then the next step is to adapt the `conf/modules.config`.

This is really where the sharing of tools on software makes it easy because we can actually share the same software but we can specify which command or parameter should be set for our specific pipeline.

So for `markdup` here, it says exactly which parameters should be used, and there a couple of different ones that I have highlighted here [7:50](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=470).

The arcs one to give you an example, and for some modules you will see have an arcs2 line. The next step (if it should only be done once), has to include this `conf/modules.config` in the nextflow.config and then at last `INCLUDE` the module into the sub-workflow or into a workflow that looks like this [8:24](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=504).

[8:24](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=504) In pink you can see the include statement that you include this module `bwa_mem` from the path where it lives, so ‘../../modules/nf-core/software/bwa/mem/main’ and so on, and then you add these options `params.bwa_mem_options` that were previously specified in the config in the `modules.config`.

So this is how you can parse them down.

Then in the sub work-flow that we have here [8:48](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=528), we have three scopes: `take`, `main`,and `emit`. `take` specifies the input data as channels, so here we have reads and indexes for my little test workflow.

`main` is where the modules come to work, so we have the `BWA_MEM` module that is included and `SAMTOOLS_SORT` modules that is also included. We can run `BWA_MEM` and take the output and run it directly in `SAMTOOLS_SORT`, and last but not least, the sub-workflows can emit named outputs to then do other things and other sub-workflows or modules.

[9:30](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=570) So here we can now name our output from `SAMTOOLS_SORT` as `sorted.bam` and then we can access this `sorted_bam` directly to do something with it.

[9:46](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=586) So last but not least, we have to combine the sub-workflows to workflows or sub-workflows of modules to workflows. So in the header, we have the two `include` statements for the workflows or the sub-workflows and the `fastqc`, i.e. once from the modules and the nf-core software, and one from the sub-workflows local.

Once again we have to add the parameters to parse them down and I have highlighted that in yellow here [10:11](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=613) for you to be able to track this.

For the FASTQC module, we just need to specify the parameters for fastqc, so modules [´fastqc´], but for the sub-workflow one needs to specify all options for all the tools in there.

[10:34](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=634) But for the sub-workflow, you need to specify all options for the tools. So here we have `bwa_mem` or `samtools_sort`, and this is exactly now this field that was originally specified in `conf/modules.config`.

[10:47](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=645) OK, then in the workflow and in nf here I have just run the module FASTQC on my input reads that I got and then the sub-workflow on the reads and the index, I get the sorted.bams output, and then I can do some more steps with that.

[11:04](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=664) So I want to show you what it looks like to add the workflows to the `main.nf` and here I’ve taken the `viralrecon` one that I mentioned earlier.

So here you see that we have three different workflows that are actually possible to do. They are apparently dependent on the input data and the whole main.nf becomes really lean. So there’s still a bit of header in here, but overall the entire main nf is less than 100 lines currently.

This makes it really easy to track which workflow is run for your input data, and you only need to look at the sra download workflow to really see what’s happening.

[11:44](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=704) OK, then there are two more important things that I would like to mention.

The first is that you always need to adapt the MULTIQC module to customise it for your tools.

This is similar to the DSL1 version that you can collect all the metrics in the `WORKFLOW` script, and then parse it down to your `MULTIQC` module as seen on the right hand side for all the different inputs. So this is one that you need to create a local module for.

You can then collect this data from the `FASTQC` module for example here.

[12:18](https://youtu.be/tWvou0xj9wA?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=738) And last but not least, you need to collect all the software versions in your workflow.

Each module must emit its software version because we need to track all the tool versions.

You can collect them all by creating an empty channel and then mixing these versions in.

From the sub-workflows, you can propagate them by these named emit versions.

So for the `bwa_mem_version`, you can get the version from the module and access it again in your workflow script as `TEST_SUBWORKFLOW.out.bwa_mem_version`, and then run your local module in the workflow.

### DSL2 - Adding modules to nf-core

Kevin Menden ([@KevinMenden](http://github.com/KevinMenden/))
will follow with a presentation on: _**Adding modules to nf-core**_

The second talk will go to the details on how to add new modules to the nf-core/modules repository and provide module tests.
