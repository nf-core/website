---
title: 'Bytesize 20: Nextflow Tower'
subtitle: Evan Floden - Seqera Labs, Spain
type: talk
start_date: '2021-09-21'
start_time: '13:00+02:00'
end_date: '2021-09-21'
end_time: '13:30+02:00'
youtube_embed: https://youtu.be/zS_hbXQmHbI
location_url:
  - https://youtu.be/zS_hbXQmHbI
  - https://zoom.us/j/94352451216
---

# nf-core/bytesize

Join us for an episode of our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize 20: Nextflow Tower

This week, Evan Floden ([@evanfloden](http://github.com/evanfloden/)) will provide some updates on Nextflow Tower (<https://tower.nf>), including the new Community Showcase with nf-core pipelines, a first look at the datasets feature, the upcoming CLI, and some automations for triggering the execution of pipelines.

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:41](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=41) What I’m going to talk about today is a little bit about [Nextflow Tower](https://landing.tower.nf/), the new features that are coming out, and some of the things we are excited about. If you are interested in a general overview, get in touch with us. We can set up a demo and also go into specific use cases.

[1:14](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=74) So here is what I will cover during my presentation. I’ll first talk about what Nextflow Tower is and discuss what it is designed to do, then I’ll cover some of the new things we have such as community workspaces, automating a lot of your workflows, reporting and outputs, and finally introducing the new Tower datasets functionality. Tower is a full-stack web application for the management of your Nextflow pipelines. It was developed after much deliberation and discussion with Nextflow users who required a couple of things. They needed a user interface to interact with a Nextflow pipeline, API actions, and a database of the history of their executions. The other thing that was really important was the configuration of environments, for example being able to set up an AWS batch to automate a process and enable users to interact the way they do with Nextflow.

[2:32](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=152) So Nextflow Tower itself is a full-stack web application. It can be deployed in your own environment or you can install it on premise as well. The application itself is made up of a couple of pieces but the main thing from the user perspective is that you can interact with it in many different ways. It has the same philosophy as Nextflow in that it can be deployed in different environments, but it can also be used on different computing platforms. So now let’s step into some of the new things we have coming out.

[3:12](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=192) The first of these is the community showcase. We really wanted to make it easy for people who are using Tower to have quick access to some of the key features around it including pipelines which have been verified. I’ll give you a quick showcase of what this looks like inside Tower.

[3:31](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=211) So anyone who signs up gets automatically added to this community showcase. This showcase is a workspace which can contain pipelines themselves that have been reconfigured. You see here that there are a couple of nf-core pipelines here, and we are adding to these all the time. All the compatible nf-core pipelines will be here over the next few weeks, allowing people to run them. This workspace is connected to a compute environment on AWS, so you can run in that environment, which is a nice and easy way to get started with these workflows themselves. There are some criteria for putting these workflows into the workspace and one of them is to ensure that there is a valid Nextflow schema file.

[4:17](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=257) So for example, if I click nf-core/rnaseq here, the fact that it has this Nextflow schema means that all the inputs can be rendered in this way. This means that it makes it easy for a user to go down and select a particular option that they want and then kick off and launch that workflow.

[4:33](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=273) The idea of this is to provide a set of validated pipelines that conform to a particular way of running. Things like `--profile test` need to be valid there, and a couple of things around how much resources these things take because we want something that would be sensible to run in a way like this.

[4:56](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=296) So now you can go here for example, and jump into that run. Since this is a shared space, I can also go look at other workflows that people have launched and see the different tasks, processes, memory etc. I’m not going to go too much into that today though and restrict myself to talking about the new stuff instead. So that was the introduction to the community workspace that we have. Another thing that I wanted to show you was the automation of this. We saw last week ([Bytesize #19](https://nf-co.re/events/2021/bytesize-19-aws-megatests)) that there was a Github action that could trigger the execution of a pipeline and that was kind of the end point there. But we’ve started to see examples of more sophisticated cases of this and they’re all visible under the actions pane.

[5:57](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=357) So you see here that we can create an action and those can have a couple of triggers. One of them could be for example a Github webhook. When I commit to this repository, this execution will fire off the pipeline - in this case a GATK pipeline. This is similar to what we saw last week ([Bytesize #19](https://nf-co.re/events/2021/bytesize-19-aws-megatests)) where we saw triggering the execution of a pipeline in some sort of a testing mechanism. What has also become really useful is the creation of a webhook here; a customised end point of an application which will allow you to trigger the pipeline if one hits the end point with some given parameters. So I can create a pipeline here by clicking on the `New Action` button.

[6:34](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=394) And then set up an action, choose a computing environment and have it all preconfigured. Then as soon as I hit the end point, the pipeline will start off. Now, this has been around for a little bit, what we are starting to see now though, is the ability to have those triggers executed based on different events. So a common one that we see with users is for a lambda function, which is able to detect a .csv file entering into an S3 bucket that will then trigger the execution of this pipeline. This is the kind of whole automation you can imagine: data coming off the sequencer, going through a pre-defined pipeline and then being notified of that. We’ve been seeing a lot of success with that as well. That was really just the first point of that. We also wanted to make it a lot easier for you to have any kind of automation here and also be able to automate different parts of the system, maybe not just launching pipelines, but other aspects of that.

[7:35](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=455) So we released the full open API for Nextflow Tower itself as the application grew to include multiple users and multiple organisations etc. We wanted to expand the API to make it fully visible to users. This allows you to then interact with it in many different ways but also allows us to build an SDK on top of the API, which allows us to build the CLI. So this is something that we have just released recently. It is still an alpha release, but I’m going to try and show you how it works.

[8:13](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=493) So you have an interaction very much similar to what you would use if you have used `nextflow run` with Tower. You essentially specify a workspace ID, specify a token that you’re going to be using, and then that submits into that workspace that you’re going to use. There are a few more things here about how you can set this up, and there’s a tutorial as well. For the purposes of this talk, I am going to give you a short demo of what it looks like.

[8:39](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=519) So if I was here in the command line interface, I could then interact with my Tower server. Maybe I’d like to list the credentials that are available to me in Tower. If I list them here, I can see I have some AWS credentials set up. When I set them up, they’ve got some ID associated with them, and then very much like `nextflow run`, we wanted to have `/towr launch`. This is for those who want to interact with the command line. I can use this to launch any Nextflow repository. If it’s set to run, I can check in Tower in the workspace to see if it has triggered the execution of the pipeline. So this is another kind of quick way to do this. We also see cases where people want to automate the creation of computer environments; one for every user or parameterise that somehow. It also allows you to do that with the CLI.

[10:08](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=608) OK, next thing is small but maybe quite useful for people who might not have outputs that are straightforward to visualise. Those could be in an S3 bucket for example. What we have come up with is a way to manage outputs, essentially to view and share results from a Nextflow pipeline. So what you can do here is define any given output (a classic example is MultiQC) in the workflow definition. Once the pipeline run is complete, you can then visualise those reports in PDF, HTML etc. formats. This keeps a kind of centralised place for those results and it allows visualisation of the different processes during the run.

[11:38](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=698) The final thing that I’m going to talk about, which is a little related to what I just said, is our first attempt to venture into data management. So I am pleased to announce the Nextflow datasets functionality. So what does this do? This provides a means to define input datasets that can then be run in a pipeline. It sounds trivial, but it is something that allows you to keep track of everything that has gone into a workflow. These can also be updated and we can keep versions of them that make it much easier for a user to be able to launch a workflow based on a given input.

[12:34](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=754) So in the datasets functionality you see an example of here, I can create a new dataset - let us call it rnaseq-example-samplesheet, and set a description or drag and drop a sample sheet as a .csv file. What I added here is the nf-core megatest. In this case, it has got a row as the header so I can specify that etc. I could reupload a new version of this and it would have a version update. The really nice thing is that this follows a scheme. So if I create this now, you can see I have rnaseq-example-samplesheet that I created earlier. Now imagine if this could be automated. As a launchpad user i.e. as someone who’d like to kick off a pipeline, I can have predefined pipelines here now, but when I go to run this input, I can simply click in the input box to see a dropdown of the different inputs. So I can select the workflow that I have and then launch the pipeline. It might seem kind of trivial, but you can imagine setting this up with outputs and then from the launch of a pipeline, you would be able to see the inputs, outputs, and be able to reuse and trace these.

[14:29](https://youtu.be/zS_hbXQmHbI?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=869) So that was it. If you’d like a more detailed demo or dive into the stuff, please reach out to us.

</details>
