---
title: 'Bytesize 2: How nf-core configs work'
subtitle: Maxime Garcia - SciLifeLab / Karolinska Institutet, Sweden
type: talk
start_date: '2021-02-09'
start_time: '13:00 CET'
end_date: '2021-02-09'
end_time: '13:30 CET'
youtube_embed: https://youtu.be/cXBYusdjrc0
location_url:
  - https://doi.org/10.6084/m9.figshare.14160347.v1
  - https://youtu.be/cXBYusdjrc0
  - https://www.bilibili.com/video/BV1M54y1a7Uy
---

# nf-core/bytesize

Join us for an episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 2: How nf-core configs work

This week, Maxime Garcia ([@MaxUlysse](http://github.com/MaxUlysse/)) will present: _**How nf-core configs work.**_ This talk will cover:

- Making your own Nextflow config file
- Changing resource requirements for a process
- Using config profiles
- How the nf-core `--max_memory` options work
- Adding a new institutional config profile

The talk will be presented on Zoom and live-streamed on YouTube:

- YouTube: <https://youtu.be/cXBYusdjrc0>

<details markdown="1"><summary>Video transcription</summary>
**Note: This text has been edited to make it more suitable for readers.**

OK, so hello everyone, Maxime here.
I am working for the Swedish Childhood Tumor Biobank which is located at Karolinska Institutet and I am sitting half-time at the National Genomics Infrastructure at SciLifeLab.
Today for the second nf-core/bytesize, I am going to talk about nf-core/configs and how they work.

What I like about Nextflow or in fact what I actually love about Nextflow, is its portability, shareability and reusability.
As a challenge I am going to run [nf-core/eager](https://nf-co.re/eager) - a pipeline I am not very familiar with - on the provided test data.
I will be assuming I am running it on a new machine with just Docker installed and I am first going to specify everything on the command-line without using any config file or any profile.
This is not something you should normally attempt.

[[1:04](https://www.youtube.com/watch?v=cXBYusdjrc0&t=1m4s)]
It is actually quite simple.
I install Nextflow with the first command.
The second command will help me download the data that we want and with the third command I'm actually running the latest release of the [nf-core/eager](https://nf-co.re/eager) pipeline.
The last command specifies all the necessary input for the pipeline:

- The container engine I want to use, which is Docker.
- The specific container I want to use is specified with a tag.
- Some resources, for example `max_cpus`.
- The path to the reference genome file, e.g fasta file.

However, we can improve this if we use config files, which is the whole point of this talk:
How do nf-core configs work?

Reading the [documentation](https://www.nextflow.io/docs/latest/config.html) of Nextflow config files, it says you can put all the parameters and all the properties that we need in the pipeline.
This simplifies the command a lot.

[[2:18](https://www.youtube.com/watch?v=cXBYusdjrc0&t=2m18s)]
So I have made a config file which I call `my_computer.config` in which I have specified how Docker should work and the resources that are available on my computer.
I do not need to specify the actual container for [nf-core/eager](https://nf-co.re/eager) because it is already specified in the `nextflow.config` file, which is provided by the pipeline.
I could also have specified the genome and the input file as well, but as I am planning to run this specific command only once, I prefer to give that on the side.

[[3:14](https://www.youtube.com/watch?v=cXBYusdjrc0&t=3m14s)]
Another benefit of Nextflow is that I can also use profiles to do the same thing.
Once again, reading the [documentation](https://www.nextflow.io/docs/latest/config.html) is the recommended way to learn how profiles work, but basically profiles are like aliases for configs.
In this new command I am using 2 profiles at once.
I am using the `test.tsv` profile and the Docker profile.

As you might have guessed, the Docker profile provides all the information related to Docker.
The test profile specifies resources for a small computer and also provides information for the input file and the reference genome.
On the contrary to the previous example, the input file and the fasta file is specified in the test profile since these will be used very frequently.

One can easily realise the benefit with profiles in Nextflow with an example.
If I want to run for example, the same command but with Singularity instead of Docker, I just need to specify the Singularity profile instead of the Docker profile.
For nf-core pipelines, the Singularity profile is available by default.
I just have to change the command to use this profile and _voilà_, it is working well!
This is why I like Nextflow and nf-core - it is easy.

So far, the examples have been designed to run on my computer. However people usually have large datasets and it is not possible to run the pipeline on a single computer. So, assuming I want to do this on my institutional server/cluster/HPC instead, we need to ask ourselves some questions:

- Which container/virtual environment engine is available to us?
- What are the available resources?
- Which scheduler or executor are we using?
- Where are the reference files?
- Where are the input files?

[[5:33](https://www.youtube.com/watch?v=cXBYusdjrc0&t=5m33s)]
If we have the answers to all of these questions, we can put it in a config file.
So for example here on this fictional server, my config will contain everything that I need to run my project.
In the Singularity scope I have defined where the Singularity containers will be located.
I have enabled Singularity because of course I want Nextflow to work with Singularity.
I have also specified some specific options for Singularity to run, so that I know that it mounts the proper folder.

In the process scope I've specified the singularity module to be loaded every time I run a job.
I've also specified the slurm executor together with some specific cluster options on my cluster to say in this case a project ID that should be used.
Finally, I specify what resources are available on my cluster.

This is now just a config file, but in order to enable other people to use this it would be good to make it into a profile.
According to the [documentation](https://github.com/nf-core/configs#adding-a-new-config) that we have on Github for nf-core configs, I should start by forking the [nf-core config repository](https://github.com/nf-core/configs).
Then I copy the config file that I have already created and I put it here in the folder `conf/my_hpc.config`.

[[7:43](https://www.youtube.com/watch?v=cXBYusdjrc0&t=7m43s)]
In the `nfcore_custom.config` file, I add a specific line that tells Nextflow to look for this specific config when I use this profile.
Of course do not forget to also update the documentation and the CI tests.
This is ideally done on a specific branch on my personal fork of this repository so that you can easily create a pull request on the [nf-core config repository](https://github.com/nf-core/configs).
Other people in the community will then have a look at it and have some comments, and eventually approve it and merge it when satisfied.

[[8:35](https://www.youtube.com/watch?v=cXBYusdjrc0&t=8m35s)]
Here are some tips:
All nf-core pipelines are designed for usage on a typical HPC, with reasonable default resources for each process.
It will actually look more or less like this for every pipeline In the `conf/base.config`.
We will have specific resources that are defined for CPU, memory and runtime in the process scope.
These are just the default ones, they are usually overridden for processes with specific labels.
There are usually several labels and we try to make them as broad as possible to work on a typical cluster but of course you might need to adjust that for your own cluster.
Finally, in the `nextflow.config` of each nf-core pipeline, the maximum resources are specified.

Another tip is that the 'max resources' is just a threshold not to go over, so changing it will change the limit but it will not change the resources that the pipeline will start with.
If you want to change the base resources you must look at the CPU, memory, and time properties in the process scope.

[[9:57](https://www.youtube.com/watch?v=cXBYusdjrc0&t=9m57s)]
Here is an example of how to change that.
Within the process scope, you can use different process selectors.
Specifically the ‘withName’ and ‘withLabel’ selectors are useful if you want to change the properties for one process or for multiple processes that share the same type.

[[10:23](https://www.youtube.com/watch?v=cXBYusdjrc0&t=10m23s)]
You can also include a config file within a profile, it could be quite useful.
You can also test your profile online if you made your PR (pull request) but if it is not yet merged, you can specify the config_custom_base pointing to your own fork of the nf-core config repository.

Some final messages to end this talk:

- Read the docs. Everything is in the docs.
- Try things out, and do not hesitate to ask questions.
- Stay tuned for more Bytesize talks.
- Get involved: Join the nf-core Github organisation, follow us on Twitter and on Youtube.

[[11:33](https://www.youtube.com/watch?v=cXBYusdjrc0&t=11m33s)]
I would like to thank everyone at my institute, sponsors and all the institutes I am collaborating with.
Here are all the institutes that are collaborating with us on this nf-core project.
And here are all the Github collaborators within nf-core.

[[11:59](https://www.youtube.com/watch?v=cXBYusdjrc0&t=11m59s)]
Finally, some concluding important links including a link to these slides.
If you have any questions, now is the time.

</details>
