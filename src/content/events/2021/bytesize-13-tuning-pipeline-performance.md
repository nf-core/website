---
title: 'Bytesize 13: Tuning pipeline performance'
subtitle: Gisela Gabernet - QBiC Tübingen, Germany
type: talk
startDate: '2021-05-18'
startTime: '13:00+02:00'
endDate: '2021-05-18'
endTime: '13:30+02:00'
youtubeEmbed: https://youtu.be/Qw1gLpYtMec
locationURL:
  - https://www.bilibili.com/video/BV1z64y1k7a3
  - https://youtu.be/Qw1gLpYtMec
  - https://doi.org/10.6084/m9.figshare.14680260.v1
---

# nf-core/bytesize

Join us for an episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation.
Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 13: Tuning pipeline performance

This week, Gisela Gabernet ([@ggabernet](http://github.com/ggabernet/)) will present: _**Tuning pipeline performance**_.

This will cover:

- In-depth analysis of computational resource usage after a pipeline run
- How to take that back into a custom config to make it more efficient

The talk will be presented on Zoom and live-streamed on YouTube:

- YouTube: <https://youtu.be/Qw1gLpYtMec>
- Bilibili: <https://www.bilibili.com/video/BV1z64y1k7a3>

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:47](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=47) I will talk about tuning pipeline performance and specify the resource requirements for your pipelines. So, here’s the background, maybe you are a pipeline developer of a Nextflow pipeline and you’d like to incorporate your pipeline into nf-core or maybe you’ve found an nf-core pipeline that you would like to use to work with your own data.

[1:13](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=73) Then you are ready to launch your pipeline to your local cluster or commercial cloud infrastructure with potentially high amounts of input data. The nice thing about nf-core pipelines is that we have a predefined set of resource requirements that should allow you to more or less effortlessly launch your pipelines across all these infrastructures. So I’m going to cover how we set this up and tweak it in case you need to.

[2:01](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=121) When writing a Nextflow pipeline, for each process of the pipeline, there is a default possibility for specifying resource requests. This is for example, the number of CPUs that this process will use, how much memory this process will need, and for how long that job can run. It is important to efficiently request those resources in the cluster or cloud scheduler. If you request fewer resources, and your job exceeds what you’ve requested, it will be killed. On the other hand, if you request too much, you end up using your cluster inefficiently, which can result in unnecessary costs. So it is important to find a balance.

[3:11](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=191) There is no magic formula to find the best configuration, but let’s have a look at something that could help. So start small and then go big. But here’s another very nice tool, which is the Nextflow execution report that can help you define your resource requests. So when running an nf-core pipeline, there’s a subfolder called `pipeline_info` under the pipeline results folder, where you can find the Nextflow execution report. We’re now going to have a look at a sample execution report for `nf-core/rnaseq`.

[4:03](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=243) So in this Nextflow report for `nf-core/rnaseq`, you can see that it’s run all the processes successfully.

[4:19](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=259) An interesting section of this report is the report on resource usage.

[4:22](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=262) Here you see that for each of the processes of the pipeline. It lists how many CPUs were used for each process.

[4:36](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=276) It also shows the statistics across the different samples that were running in the pipeline.

[4:46](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=286) You can see the raw usage and there’s another tab that allows you to see what percentage of the allocated resources were actually used.

[4:57](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=297) In this case, you’d like to optimise that so that you use most of the allocated resources, while also leaving a margin for larger samples that can take longer than usual.

[5:16](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=316) We also have a section for memory requests. Here you can see the physical RAM that each of the processors used, the virtual RAM, and the percentage of the allocated RAM that was actually used by the processors.

[5:36](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=336) Finally, there’s also a section on the duration for each process.

[5:48](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=348) One recommendation would be to use some defaults that seem sensible and then have a look at the resource usage report.

[6:10](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=370) When you’re developing a pipeline and writing the code, this can be specified as a part of the Nextflow configuration, so that would be a part of `nextflow.config` for most Nextflow pipelines. For nf-core pipelines, we have this resource configuration in a separate file, so that’s easier and cleaner to read. That one is called `base.config` and can be found in the subdirectory confs of the pipeline. This base.config is imported inside the main nextflow.config, so inside this config file there’s the process scope where the defaults of CPU, memory, and time requested for each process is defined. In this case, the default would be one CPU, 6 GB of memory, and four hours of time. This is not going to be enough for all the processes in the pipeline. Another thing that we can see from this default configuration is that we have these resources multiplied by `task.attempt`. Whenever a process is retried, the number of requested resources will be doubled. These can be defined with the exit statements `errorStrategy`, `maxRetries`, `maxErrors`. We can also see if a process should be retried or finished depending on the exit status. So in this case for nf-core/rnaseq, if the jobs end with one of these exit statuses, that process will either be retried after doubling the initial resource requests or will be finished with this `maxRetries` option. You can also define how many times the process should be retried. In this case, it’s just one, so it will be retried once. If it ends with one of these exit status codes, it will be retried with double the number of requested resources. Finally, another thing that we see here in the default process definition is this `check_max` function. This ensures that the resources you request for each of the processes don’t exceed the available maximum computing resources. We will later see how to define these, but in a nutshell, if you’d like to launch a pipeline within your institutional cluster, you have a maximum amount of memory and number of CPUs that each node has. So this ensures that this request on exit does not exceed those maximum resources.

[9:28](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=568) It could be the default that is there for each of the processes, but these would likely be insufficient for all the processes in the pipeline. So it should be defined if some processes need more resources. For nf-core pipelines, we define this in the base.config that specifies different types of processes. So for example process_low that will need more resources as a default but still low resources or process_medium that will increase to a minimum resource definition here. This is done via labels and this code will also be a part of the nf-core template, so you don’t need to define this. You should just try to use those labels as they are defined in the base.config when defining your own processes. So for example, for the `FASTQC` process, this one contains a process_medium label already, meaning that for this process, these resource requests will be used. You can also define these labels for new processes as part of nf-core pipelines. Then the resource request definitions are already pretty fine.

[10:58](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=658) However for some of the processes, you might need to define some defaults that don’t match any of the labels here. In this case you can also define the base resources request for that process `withNAME` statement, and there you provide the name of the process directly instead of a label.

[11:26](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=686) In the case of DSL2, the name of the process should still be the name that defines part of the pipeline code and not also include the name of the workflow or sub-workflow where this module is present. So that’s how a specific request can be defined.

[11:56](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=711) So that in principle is taken care of for nf-core pipelines. However, sometimes you might want to change those defaults that are defined for a specific pipeline. So when you’re using an nf-core pipeline, one of the first things you need to define resource-wise is what the maximum resources you have are, be it in your local computer or on the cluster that you are using. This can be done directly via pipeline parameters - all nf-core pipelines have the parameters `--max_memory`, `max_cpus`, and `--max_time` - and you can directly specify the maximum resources there. But if you’re using nf-core pipelines often within your cluster, it is a good idea to define an institutional config profile (see [Bytesize#10](https://nf-co.re/events/2021/bytesize-10-institutional-profiles) as part of the nf-core configs repository. The max resources would already be defined as part of this institutional profile.

[13:24](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=804) What if you would like to define not only the maximum resources but also tweak the resource request for a specific process because you have especially large samples or small samples, you can do that with a custom nextflow.config. It’s also possible to provide nextflow.config with the -c parameter, and inside this custom config, you can define using the same syntax as before. Special resource requests for the process that you’d like to change for example in the case of MarkDuplicates, the process requests a change to 20 GB and CPUs to 20. So that’s in case you’re a pipeline user and you’d like to change requests for a specific process. If you have that often, you would need to change the defaults of a specific pipeline when running it for your data. You can also contact us on Slack so that we take care of these defaults and tweak them so that most users don’t face those issues. Or even open a pull request yourself to the base.config to change these defaults.

[15:09](https://youtu.be/Qw1gLpYtMec?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=909) So that’s what I wanted to cover today. Reach out on Slack if you have any questions.

</details>
