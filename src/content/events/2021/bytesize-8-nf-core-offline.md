---
title: 'Bytesize 8: Running pipelines offline'
subtitle: Maxime Garcia - SciLifeLab / Karolinska Institutet, Sweden
type: talk
start_date: '2021-04-13'
start_time: '13:00 CEST'
end_date: '2021-04-13'
end_time: '13:30 CEST'
youtube_embed: https://youtu.be/N1rRr4J0Lps
location_url:
  - https://youtu.be/N1rRr4J0Lps
  - https://www.bilibili.com/video/BV1AZ4y1w7hZ/
  - https://doi.org/10.6084/m9.figshare.14413046.v1
---

# nf-core/bytesize

Join us for an episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation.
Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 8: Running pipelines offline

This week, Maxime Garcia ([@maxulysse](http://github.com/maxulysse/)) will present: _**Running pipelines offline.**_
This will cover:

- Using nf-core `tools` to download a pipeline
- Download Singularity images, with and without nf-core `tools`
- Download DSL1 and DSL2 pipelines

The talk will be presented on Zoom and live-streamed on YouTube:

- YouTube: <https://youtu.be/N1rRr4J0Lps>
- Slides: <https://doi.org/10.6084/m9.figshare.14413046.v1>

<details markdown="1"><summary>Video transcription</summary>
:::note
This transcript has been modified to make it reader-friendly
:::

[0:32](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=32) Thank you, I will share my screen.

[0:40](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=40) So hi, Maxime here. Today I’m going to talk to you about how to run a pipeline offline. I’m going to try to run two different pipelines.

[0:56](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=56) So the plan is rather simple. First, we need to download the DSL1 pipeline, so we do that on `sarek`, the main container and an extra container because it’s `sarek`.

Then I’m going to download the DSL2 pipeline using `nf-core/tools` to get the container for all the modules, just in one command.

Then we need reference data, so I’m going to download AWS iGenomes.

Then I will download Nextflow as well.

Then I’m going to transfer all that to an offline machine and then run the pipeline.

I’ve done the downloading already, so I’m not doing it here. This is not a live demonstration.

[1:45](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=105) So my set-up is loosely based on the system I am using.

In Sweden, we have this server, which is completely offline and we can send files to the server via a shared folder that we can connect to with `sftp`.

That’s what I will be doing now.

[2:11](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=132) So first on the machine that can actually communicate with this offline machine (I am connecting to this machine) with the `ssh` command.

Then I’m going to install the latest version of nf-core/tools.

I’m using the `dev` branch, so that I have access to everything.

You can see I have installed version 1.14.dev0.

[2:43](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=163) So now I’m downloading the DSL1 pipeline, so `sarek`.

I’m specifying a version, 2.7, so that’s the latest version, and I add the `-s` flag, which means that I also want nf-core/tools to download the `sarek` container.

[3:13](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=193) So here we can see that `nf-core/tools` has found three containers, but it’s downloading just one.

Because the extra containers are specified in the `config` file, they are not recognised by `nf-core/tools` but this will be solved in DSL2.

[3:32](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=212) So here I am using the simple Singularity command, so `singularity pull` to download the extra Singularity container.

So that’s fairly simple.

With this command, I can also download the offline Singularity container for the pipeline.

I could show that, but it’s the same command. It’s easy to split it and download the main one with `nf-core/tools` and download the extra one with `singularity pull`.

[4:14](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=254) Downloading a DSL2 pipeline is simple too.

This slide shows how to download the `nf-core/rnaseq` pipeline.

We need to specify the release, which is 3.0.

I specify this with the `-s` flag that I want to download the Singularity image.

I have access to a very big machine.

I use `-p10`, which means that I will download 10 Singularity images at once (in parallel).

10 is a good number for an example although this machine would allow me to do more.

If you do that, you will see a progress bar, but this process takes a bit of time, so I’m not going to show that here.

[5:19](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=319) I also use `--compress none` because I don’t want the final folder to be compressed because with that many containers (29 here), it will take a long time to compress them all.

Besides, I won’t really save much space.

I think it’s good to compress if you just have the script, but if you have you have dozens of containers, I don’t think it’s worth the effort.

[6:12](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=372) Then let’s try and download the rest.

So first I want to download Nextflow.

If you look at the Nextflow release, you have the Nextflow download URL with the name of the release, and then `-all`.

If you don’t know this one directly, you don’t need internet access for Nextflow, which is very useful in our case.

In my case, I have access to Nextflow through the cluster I work with, but in this current fictional case I don’t, so I’m going to download Nextflow myself.

In real life I also have access to the AWS iGenomes set up, but again, let’s assume it’s not set up.

So I need to download it first so that I have it on my machine as well.

So I’m using the simple `aws s3` command to download everything from the s3 folder.

[7:31](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=451) Here I’m downloading the GRCh folder for sarek and nf-core/rnaseq.

We have two different versions because for `sarek`, I’m using the GATK bundle.

[7:50](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=470) And then I connect to my offline server using `sftp`, and then the simple `sftp` command.

So I can use `-l -ls` to look at all the files I have in local, and then I use `put` to move the file from the local server to the distant server.

I move the whole folder using `put -r`.

That transfers all the files.

Then you need to connect to your offline server.

Now this is not strictly offline, but what I mean is that it’s a server with an open port and through which you can’t connect to anything else.

On this offline server, I move my file around. For example I put the Nextflow -all in bin, I make a symbolic link for that and I also put my pipeline in my own directory, so I know where they are.

But that’s more or less easy to do.

[9:26](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=566) Then it’s fairly simple to run the pipeline.

So here I’m going to show what you can do if you don’t really use a profile.

So I’m specifying the Singularity cache directory with `$NXF_SINGULARITY_CACHEDIR`.

So you can do that, but I think it’s better if you specify the cache in the config file.

[10:27](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=627) So the first example will be how to run the pipeline without specifying much in the profile, and the second example is running the rnaseq pipeline.

I suppose that everything has been specified in the offline profile.

[10:46](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=646) That’s all actually.

So here are some tips.

If you have a different set up, for example if you have access to a head node that is online and the others are offline, you can set up `$NXF_SINGULARITY_CACHEDIR` to specify a central cache for container downloads.

The latest version of nf-core/tools is quite smart, and won’t download the same container more than once.

Then use parallel for faster downloads, and `--compress none` because it’s much simpler not to compress big folders.

[12:05](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=725) Apart from that, my main tips are read the docs, try things out, and don’t hesitate to ask questions.

> **Conversation after comment from Phil**

[13:44](https://youtu.be/N1rRr4J0Lps?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=824) Just at the end of your slides Maxime, you mentioned that you needed to set the config base parameter, is that right?

Yes.

OK, that used to be the case, but that’s not true anymore.

Ah yes, I used that because it was in the `README`.

We need to update that then, sorry.

So if you use `nf-core/download`, it will create the directory with three sub-directories i.e. the workflow, the singularity images if you use `-s`, and the config.

The pipeline files are actually edited so that they know where everything is.

So you don’t need to do anything.

You just do `nextflow run` in the workflow directory and it will know where the config profiles are automatically.

OK, that’s good. No I noticed that there was this kind of architecture, but read the `README` and I’ve been doing that so that’s why I thought that that was how it still worked.

I’m glad to know that it’s now possible.

Sorry, we just need to update the docs in that case.

Ah no problem, now we know.

Thanks.

</details>
