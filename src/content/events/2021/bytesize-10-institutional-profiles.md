---
title: 'Bytesize 10: Making a new institutional profile'
subtitle: James A. Fellows Yates - LMU Munich / MPI-EVA
type: talk
startDate: '2021-04-27'
startTime: '13:00+02:00'
endDate: '2021-04-27'
endTime: '13:30+02:00'
youtubeEmbed: https://youtu.be/Ym1s6sKGzkw
locations:
  - name: Online
    links:
      - https://www.bilibili.com/video/BV1BA411V78U
      - https://youtu.be/Ym1s6sKGzkw
      - https://doi.org/10.6084/m9.figshare.14541291.v1
---

This week, James Fellows Yates ([@jfy133](http://github.com/jfy133/)) will present: _**Making a new institutional profile**_
This will cover:

- what information to gather before making a institutional profile
- step-by-step construction of an example institutional profile
- how to test a draft institutional profile

The talk will be presented on Zoom and live-streamed on YouTube:

- YouTube: [https://youtu.be/Ym1s6sKGzkw](https://youtu.be/Ym1s6sKGzkw)
- BiliBili: [https://www.bilibili.com/video/BV1BA411V78U](https://www.bilibili.com/video/BV1BA411V78U)

You can see the slides on HackMD: <https://hackmd.io/@jfy133/Hyc_WvXU_> (shown below)

<div class="ratio ratio-16x9 border shadow">
  <iframe  src="https://hackmd.io/@jfy133/Hyc_WvXU_" allowfullscreen></iframe>
</div>

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:54](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=54) Thank you very much. As mentioned, the overview of today’s talk is that we will be recapping the benefits of using nf-core institutional profiles, and I will describe the information you should generally gather before writing one.

[1:14](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=74) Then I’ll go step-by-step through how to build such a profile, and finally how to test it.

[1:20](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=80) So today, I’ll be talking about global institutional profiles, which is what I will be defining during the talk. So institutional profiles are Nextflow configuration files that are used by all nf-core pipelines to work efficiently on institutional level clusters and that are stored on the `nf-core/configs` Github repository.

[1:40](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=100) There’s another thing called pipeline institutional profiles, but that is something for a future bytesize talk.

[1:47](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=107) So firstly, why would you want to use an institutional profile?

So these give you a lot of efficiency regarding Nextflow runs in terms of computing resources and time.

So making sure that you use only what you need and not exceed it and block other people in your given cluster.

[2:04](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=124) They’re very portable and contribute towards reproducibility, so that means other people can run the same command you ran, but configured to their own cluster.

[2:15](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=135) So ultimately writing this `config file` saves a lot of time. You write it once and then everyone on that cluster can benefit from it.

[2:25](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=145) So here’s a brief recap of what was covered in [bytesize#2](https://nf-co.re/events/2021/bytesize-2-configs).

A Nextflow configuration file is a simple text file containing a set of properties, so parameters and some stuff that basically Nextflow can read and understand. This can be pipeline parameters, information on locations of files, and so on.

[2:49](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=169) There are different levels of these configuation files, so whenever a Nextflow run starts, it will look in the `nextflow.config` file that may be present in your run directory.

It also checks your home directory for the `nextflow.config` file.

You can also specify custom config files with `-c` parameter, and also in nf-core pipelines there’s a whole other layer, which are nf-core profiles, which like I said before are sorted `nf-core/configs`.

This is generally where we store our institutional profiles.

[3:22](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=202) While I say that this is stored on the Github repository, that’s not to say that you can’t run these offline.

So if you have a cluster that is on offline mode, you can also use `nf-core/tools` to basically clone these config files to your cluster.

[3:40](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=220) So here an example of what exactly I mean by these profiles and how you can use them.

Let’s say I want to run `nf-core/eager` on the UPPMAX cluster in Sweden, I can simply write `nextflow run nf-core/eager` and then specify the profile as Uppmax i.e. `nextflow run nf-core/eager -profile uppmax`, which is the institutional profile, and then all the other parameters I want to use.

But as I said, these profiles can apply to all pipelines, so you can also switch at the bottom to `nf-core mag`, and the same thing will work, but everything will be optimised for that particular cluster.

[4:12](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=252) You can also chain multiple institutional profiles together if you need.

So for example for `nextflow run nf-core/eager -profile shh,sdag` where we have, dash, then the profile for the MPI in Jena, we have multiple clusters.

So we can specify within that specific cluster the specific institution that you want to use.

Something that’s very important that was mentioned in [Maxime's talk](https://nf-co.re/events/2021/bytesize-8-nf-core-offline) is that the _order_ of the profile is very important.

So anything on the furthest right will override anything previously. So any parameters set in the `sdag` profile will override any parameter also set in the `ssh` profile in this given example.

[4:57](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=297) So let’s say you want to write an institutional profile for your cluster, there are multiple things you may want to cover.

[5:05](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=305) And this is the information that you want to gather before you start writing the profile. So it will cover things like names, resource limits on scheduling systems, and also containers.

[5:15](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=315) So the first question you want to ask is what name do you want to call your profile.

This should be something that’s very recognisable and also descriptive so that users of your cluster can recognise which cluster is being referred to.

But it should not be something that is too generic either such that people might confuse it with their own cluster.

Short is also very good, abbreviations are OK, like the `shh` example from earlier.

But again, it should still be precise enough for people to understand what it is.

If you have multiple clusters for one institution, it might make more sense to have an umbrella name for the institution itself rather than have one profile for each cluster.

Within that you can have internal profiles, which I’ll give you more examples of later on.

So make sure you get the right level of hierarchy there.

[6:13](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=373) Once you decide on a name, you need to go into your cluster and get into information about the specifications of that.

For example, you should look for the resource limits of that particular cluster.

So nf-core pipelines expand on the Nextflow retry system where if a process times out because it runs out of memory or time, it makes it re-submit that job with more memory or time.

However nf-core adds to these limits to ensure that this retry system doesn’t start requesting more than is actually physically available on the system.

So for this, you need to look up the largest node that all users of your particular cluster can access, and find out what the most amount of memory, CPUs and walltimes are.

So for example, a 2TB node with 112 cores and the longest walltime it can have is seven days.

In addition, you should also look up whether scratch usage is required by your cluster (this is a space that doesn’t get backed up, where all your working directories and files go into).

[7:23](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=443) So once you’ve found that you need to look into your scheduling system, assuming your cluster has one (not all do..), and check that Nextflow supports that.

So a scheduler is basically a tool which allows the computer to distribute jobs where resources are available for you.

This can be based on priority, so you should look up if your scheduler has queues or partitions that are often based on walltime, meaning that shorter jobs will go earlier in the queue whereas longer jobs go further down.

You should find out the names of these and what specifications you may have for those, you should also look if there are any submission limits, so for example you can only submit x number of jobs at any one moment.

You should also check if there are any additional configurations that you may have when you’re writing your own custom batch scripts that you would feed to the scheduler.

So something like `module load xyz`.

[8:35](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=505) Then finally, you should look into containers.

So nf-core highly promotes the use of container systems.

This is stuff like Docker, Singularity and so on, which basically allow you to put all of the software that a pipeline needs into a single file.

Then you just download that file with all the tools set up and configured in the way it needs to run successfully.

That basically allows for much more robust reproducibility because everyone will always use exactly the same parameters and so on.

So you should check whether your cluster offers the use of these, and if so, which container engine, and then check if nf-core supports this.

Remember that nf-core has specific profiles for this stuff like Docker and Singularity, Charlie Cloud, Podman and a couple of others.

You should check if your cluster has a common cache location; this is a directory where all the images are stored in one place and that all users can access.

This prevents loads of the users having the same copies of the image.

You should also check whether you need any additional parameters to these container engines, for example bind paths and so on.

[9:35](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=575) so once you have all this information, you can get started with writing.

[9:40](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=580) We recommend starting with making a fork of the `nf-core/configs` repository.

You don’t have to do this, but it generally makes life easier when you want to basically image everything.

You can make a branch to follow GitHub best practices and just name it after the cluster name you’ve picked.

Then also open a tab with the nf-core and Nextflow documentation because you will have to refer to it quite a lot.

[10:03](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=603) Once you’ve made your fork and your branch, you need to create two new files.

One under the `configuration` directory, which is where you’re going to put in the configuration file itself, a documentation file in markdown, and docs.

You then need to add your profile name to these three files that nf-core uses to organise its stuff.

It’s quite straightforward where to put your name there, just look in and you will find where the long lists of the other profiles are.

[10:32](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=631) Now to write the profile itself. So in the configuration file, you can start by writing the `params` scope, this is basically standard Nextflow parameters, which you often specify in the command line.

[10:44](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=644) So for nf-core profiles, we generally suggest that you add to these the config profile description content and url. This is useful later on for debugging, and it also makes it easier for users to know that they’re using the right profile when they are running their pipelines.

[10:59](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=659) Then in the `max` parameters that you see in the 4th to the 6th line, you can specify the maximum and node information that you gathered earlier. So the maximum number makes for memory that’s available, number of CPUs, and so on.

[11:14](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=674) If you for example have a common directory for things like reference files, which you use in the network pipelines themselves, such as the Illumina iGenomes resource, you can also specify this here.

[11:30](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=690) Once you've specified that, you can then write the `process` scope.

So this is actually where you define generic information about your particular scheduler.

So in the example above (see Simple Example on the slide), we see the setup that some clusters can have, and you just specify using slurm and maybe add an additional security check so that you can only retry a maximum of two times to stop runaway processes before you hit the big resource nodes.

[11:55](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=715) In a more complex example (see Complex Example on the slide), you can also specify queues and do this dynamically.

So using a `groovy` expression, you can say for example if the process we’re going to request is going to run for less than two hours, we go to the short queue.

If not, for less than 24 hours it goes to the medium queue, and for everything else, it goes to the long queue.

[12:24](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=744) Then also, if, for example, your batch scripts for your scheduler that you normally write manually, you have to specify extra parameters in that header block (you can also see this in the cluster options), for example saying that in my SGE batch script I need to specify the `h_vmem` parameter.

Then I can take the Nextflow task memory information using the groovy curly brackets.

[12:55](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=775) Then you can specify more information about your particular scheduler in the executor scope.

This typically just limits queue size, which is maybe what you want to have, when you’re running eight processes for any given pipeline run.

Say you only want to submit ten processes per second, you can limit that there as well.

Again, I highly recommend checking the Nextflow documentation. There’s many more options that you can put here, and there are just a couple of examples.

[13:29](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=809) There is typically one container scope per container system.

So for example there will be one for Singularity, another for Docker, one for Podman, and all of these have their own settings.

But common ones are for example, as a non-expert user for a specific container, you can enable that with the `true` variable there, and you can specify cache directory here.

So this is if you do have that common directory where all people’s images go, and everyone can reuse the same files.

You can specify that in this scope.

[14:04](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=844) Finally, in our configuration file.

If you have multiple clusters within your particular institution, you can specify this in the `profile` scope here.

So for example in my sample institution here, I’ve got the red cluster and the blue cluster, and you can see [here](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=865) that they have different max CPU parameters.

This is how you specify that these will basically overwrite anything specified in the base institutional profile.

Again, you should add the description just specifying which cluster you are using, so that people can check and load the correct one.

Alternatively, you can do some magic, click on the link below on host names (check [slides](https://hackmd.io/@jfy133/Hyc_WvXU_)).

Use the host names for example in the UPPMAX profile, which allows you to dynamically set the parameters of different clusters based on the host-names of the machines themselves.

[15:07](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=907) Then moving on to documentation.

In the markdown file, there’s a range of things that you should include.

You should also give a brief description of where the cluster is based so that people that they are accessing the right one.

You can give a summary of the parameters that you have set, just so they know what to expect when they run the profile.

So for example, the resource limits or the queues.

[15:25](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=925) You should also write instructions for user-level configuration.

So for example if there is no common cache directory for all users on the containers, and you want to specify this at a user level, you should describe this here.

You can also describe what other sub-profiles or internal profiles of you particular institution may have.

[15:46](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=946) Once it’s all ready, we can test and submit.

[15:51](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=951) So to test the profile from your fork, you can use the special nf-core flag `--custom_config_base` and basically give this raw URL replacing your GitHub and then your branch.

Then just run this command `nextflow run nf-core/<fav_pipeline>`, and using the test profile we recommend to see if it runs.

Expect a trial and error here since every cluster has their own special quirks.

It can take a bit of time, but don’t give up.

Once you get there, you’re much more familiar with the cluster, and it’s much more straightforward.

Things you should look for when testing are things like checking to see if you see the correct description in the run header.

So this description URL and contact of your particular profile gets displayed in the nf-core run header, so this is where the ascii logo is.

If you see this, you’re on the right track.

You should also check when you’re running the pipeline to ensure that your jobs are being submitted to the queue you’re expecting.

Look to see if your jobs are being listed in your job log history, so for example in Slurm that’s with `sacct`.

You can just check that they are being displayed and so on, because if they’re not, you might still be running on the head node, which is a bit scary.

Finally, you can also check your cache directories to check if you have your containers in there as you specified.

[17:19](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=1039) So once that’s all running and you’ve tested everything, you can then submit to `nf-core/configs`, make a pull request and ask Slack for a review.

Once it’s approved, you can merge it in, and spread the word to your colleagues, so that users can get started.

They simply have to run `nextflow run`, their pipeline, `-profile` and the name of their cluster with the additional parameters you may have.

[17:49](https://youtu.be/Ym1s6sKGzkw?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e) So to conclude, if you need any help, you can always ask questions, and check the documentation.

Thank you.

</details>
