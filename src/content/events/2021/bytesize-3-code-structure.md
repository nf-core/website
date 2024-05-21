---
title: 'Bytesize 3: Pipeline code structure walkthrough'
subtitle: Gisela Gabernet - QBiC Tübingen, Germany
type: talk
startDate: '2021-02-16'
startTime: '13:00+01:00'
endDate: '2021-02-16'
endTime: '13:30+01:00'
youtubeEmbed: https://youtu.be/FFTNVbdD5pQ
locations:
  - name: Online
    links:
      - https://doi.org/10.6084/m9.figshare.14160677.v1
      - https://youtu.be/FFTNVbdD5pQ
      - https://www.bilibili.com/video/BV1gp4y1h7G8
---

This week, Gisela Gabernet ([@ggabernet](http://github.com/ggabernet/)) will present: _**Pipeline code structure walkthrough.**_

She will talk through the different files found in the nf-core template and describe what they do.

The talk will be live-streamed on YouTube:

- YouTube: <https://youtu.be/FFTNVbdD5pQ>

<details markdown="1"><summary>Video transcription</summary>
**This text has been edited to make it more suitable for reading**

Welcome everybody, I would like to do a pipeline code walkthrough today of the DSL1 pipeline template, and also provide a preview of the DSL2 nf-core pipeline template.

All the nf-core pipelines are created by using a pipeline template that contains most of the functions that the nf-core pipelines have. This makes it easy to automatically update the nf-core pipeline code with the newest developments.

To generate a new nf-core pipeline, you first need to install the nf-core tools with pip or Conda and then run the command nf-core create.

Doing this automatically creates a pipeline template, and you will first be prompted for some pipeline information: like what you would like to name your new pipeline, a short description of the pipeline, who the main authors of this pipeline are, and in less than a second, you can have a pipeline structure ready!

However, as you might see, the first pipeline template already contains six directories and 26 files that are visible. Now this together with hidden files can be around 40 files, so it can be a bit daunting to get started with developing or contributing to nf-core pipelines.

I'm going to demonstrate the function of all of the individual files inside the template, and talk about why we have them there.

We will start easy.

First of all, the most important file of the whole template is the main.nf file, and this contains the main Nextflow pipeline code. No further introduction is needed for this one.

Additionally we have other files, and different Github repositories for different projects such as the changelog/markdown file. This one contains a summary of the changes that have been introduced into the pipeline between the different releases.

We also have a code of conduct file that describes the expected behavior of people that want to contribute to the community and to the pipelines.

We also have the license file, which describes the MIT license that all of the nf-core pipelines adhere to or apply, and a README markdown file.

The README markdown file contains the easiest pipeline recommendations, such as how to get started to launch the pipeline.

What are pipelines actually? They enable the execution of different tools one after another, so it's pretty important for any pipeline to have proper software packaging possibilities. This is natively supported by Nextflow; Nextflow natively supports the execution of pipelines using many different container engines and also Conda.

We support this also of course in all of the nf-core pipelines.

So how is software packaging defined?

As part of the nf-core pipelines, we have mainly two files: the Dockerfile and the environment.yml. The environment.yml is a Conda environment file that lists all of the software dependencies of the pipeline and their versions. It can be used to just execute the pipeline by using Conda.

However due to reproducibility issues that the conda environment can sometimes have, we recommend actually running the pipeline in a containerized fashion instead. This can be done mainly with Docker containers or Singularity containers. All pipelines support this, and there are other container engines that will also be supported.

And what is the Docker container like?

So for all the pipelines, there is one Docker container available that contains all the dependencies, and I'm going to show you how this container is defined in the Docker file. We have pretty minimalistic Docker files in the nf-core pipelines, so the Docker file takes a container image that is directly provided by nf-core tools. What this does is that it copies this environment.yml file there, and creates a Conda environment inside of the container that adds this environment to the path.

So you can see that all the tools for running the pipeline in the containerized fashion are also defined in this Conda environment file. This adds a pretty strict requirement of all of the nf-core pipelines; all of the software is available in Conda, so - at any of the Conda public channels. We have two main directories here as part of the pipeline template and one directory is called /bin, and this contains executable scripts that are used as part of the pipeline and can be written in any language. Here for example, it is in Python. Adding scripts to /bin provides the additional advantage of Nextflow adding all of these scripts to the path directly. So when you want to call descriptions as part of your pipeline, you can call the script directly and Nextflow will find it inside the bin.

Additionally, we have another directory that's called `/assets`, and this directory contains the different templates that are used to send emails at the end of the pipeline saying what the status of the execution is: why the run passed or the run failed etc.

They also contain a file that is called MultiQC config and this file exists because all of the pipelines (or most of them) use the multi-qc tool to aggregate the quality control results for the pipelines and report it. This MultiQC config defines how this report should be displayed.

So that is most of it.

The assets directory also contains the pipeline logo. The pipeline logo will also be automatically created for you once you generate the pipeline with the nf-core template.

All nf-core pipelines are also required to have extensive documentation, so the docs directory contains extra documentation in addition to the main README file of the pipeline. The two most important files are the ‘usage markdown’ file and the ‘output markdown’ file. The former describes exactly what is needed to run the pipeline and provides a description of different options that the pipeline can use. The ‘output markdown’ file describes the output that we expect from all of the nf-core pipelines. So these documentations files can not only be read on the Github repository but are also automatically displayed on the nf-core website.

In the nf-core website there is a page for each individual pipeline, and you can find some tabs on display on this page containing first of all the README markdown, and also the usage docs and output docs of the pipeline. If you look closely, you will notice another tab called ‘parameter docs’, and that displays all of the parameter options that each of the nf-core pipelines take. Originally this parameter documentation was also included in the usage.markdown, but now we have an additional characteristic of all of the pipelines and that is the Nextflow schema defining all of the parameters that the pipelines take, accompanied by a short description.

This schema is then parsed as I mentioned for the nf-core website. This is not only used to display the documentation for all of the parameters on the nf-core website, but will now also be used to validate the parameters that are passed to the nf pipelines, and used to launch new pipeline runs with a parameter json file.

However nobody wants to directly modify a .json file and add the description by hand. There is an nf-core tool that helps you build the parameter schema. If you run or build a schema in the current directory where you have the pipeline, this can launch a command line client or a web interface that will help you edit the parameter schema for the pipeline and display it nicely for you to edit.

So that was the documentation. We go to the next big thing in the pipeline structure and these are the config files.

Config files can be daunting, and Maxime already introduced some of the possibilities to provide configurations with your pipeline (check Bytesize Talk #2). But I want to summarize the main files that we have here as part of the pipeline template.

The main file next to a configuration file in a pipeline, is the nextflow.config file and this includes all of the other configurations that we also find under the config directory file. The main nextflow.config file also contains core profiles that define the different parameters that can be used to launch the pipeline with managing the dependencies with Docker or with Conda. We also have different test profiles that I will describe a bit better later, and the nextflow.config files also include all of the server profiles.

Maxime described them last week - check Bytesize Talk #2; server profiles can be used to define the hardware requirements of clusters and different servers that can be used to run nf-core pipelines at different institutions. We have them specified in an independent nf-core repository that is called nf-core/configs. So in this way, all of the pipelines can import these hardware configuration files and they can be used to run any pipeline in nf-core.

Additionally in this nextflow.config, the base.config and igenomes.config files are imported, but they have been added here in a different directory just for clarity. So the base.config file contains the the default memory, CPU and time requirements for all of the pipeline processes, and it can be specified for specific processes if they need more resources or less resources.

Then we have an additional config file that is called the igenomes.config, and this one is used to specify the paths to the reference data that are needed for running most of the nf-core pipelines.

So the FASTA genome or the vba index files with igenomes can be a bit complicated; we will probably do a separate Bytesize talk to explain all the details of igenomes.

Finally, as part of this conf directory, we also have the the test profiles that are of course also imported as part of the nextflow.config.

So what do the test profiles specify?

The test profiles specify the path to the test data. The small test data can be used to quickly test the pipeline and for our integration tests. So the small test profile as I call it, points to data that is stored in a separate repository in the nf-core test data set repository, and can therefore be used by all of the pipelines.

Additionally, we now have a full test profile, and this aims to specify the path to full size data. Because of the size of this data, we do not have it available anymore in the nf-core test dataset profile. But, the idea is to use this full test profile data that is available in public data repositories such as SRA or ENA, so this data is then directly pulled from these publicly-available repositories and can be used to run the full pipeline tests. So here in the test profile is where the path to this test data and to default parameters is specified.

But how are continuous integration tests actually run with Github actions in all of the pipelines? I'm going to show you the hidden files of the pipelines, and that's where all of these Github-specific parameters and also the Github action workflows are specified.

So when we check the hidden files of the pipeline, we have a main folder that is called .github and this specifies important things such as where to interact with the Github web interface for example. It also contains issue templates to help create new issues, bug reports and feature requests for all of the nf-core pipelines. And, it contains a pull request template to help you create and go through the steps of contributing to a pipeline. As part of this .github repository, we have also a sub directory that is called workflows, and this is where all of the github actions workflows for the pipelines are defined.

So we use github actions, which are directly provided by github to trigger our continuous integration tests. The github actions workflows can be triggered at different events for example when we push to certain branches of the pipeline, when we make a pull request, or on pipeline releases.

As part of the continuous integration tests we have in several workflows. One of the workflows is called linting.yml, and it performs a markdown linting test to make sure that all of the markdown files in the pipeline follow markdown standards and nf-core linting tests. The end of call linting tests run the nf-core lint command, and this ensures that all of the nf-core pipelines follow the template standards.

Then we also run the actual pipeline tests with a small test data that I've described before. This runs a pipeline with a `-profile` test parameter. There are also two other workflows that are important for the continuous integration tests of the pipeline that are called 'push docker hub'. What these workflows do is that they build the docker container using the environment.yml file and the docker file of the pipeline and push this container to the nf-core docker hub repository. There are two different workflows, one for building the developer development container and one for building the release container of the pipelines.

Additionally, since recently we also run the full-size tests for each of the nf pipelines, and those are run on AWS (Amazon Web Services). There is also a specific Github actions workflow that is called AWS full tests, and this workflow will basically launch the pipeline full-size tests on AWS batch. One last workflow is the branch.yml workflow which ensures the branch protection restrictions that we have in the different branches of the pipeline repositories.

So that is basically it for a quick run of the DSL1 template. However in July 2020, Paolo announced Nextflow DSL2, and you might be wondering why are you even explaining to us now the DSL1 template when we are going to change things pretty soon to the DSL2 template.

The truth is a lot of features of the DSL1 template remain, so I am just going to give you a sneak preview now of the DSL2 template and you will see the main changes to the DSL1 template. To have a sneak preview of the DSL2 template, you just need to install the development version of the nf tools and check the code at the DSL to check out the branch DSL2_template and install the tools there. Then create a template with this nf-core and tools version and you will have a look at the DSL2 sneak preview that I will show you now.

So most of the template structure stays the same in DSL2, but there are some main changes I am going to describe a bit here. First of all, with DSL2 come modules, and Nextflow modules will be contained in a sub-directory called modules. There, we will have both local pipeline modules and also nf-core modules that can be shared across several nf-core pipelines because modules will specify their own containers.

Also the Docker file and .yml files will no longer be needed, and additionally in the DSL2 template, there will be a subdirectory called lib that contains different groovy functions and classes that can be imported in the main.nf file so that makes the main.nf file much cleaner. So that most of the boilerplate that is now in the main.nf file will be a part of the lib sub-directory. Additionally we will have a modules.config that defines the different parameters that will be needed to run each of the modules of the pipeline. But just be warned; this template is unreleased yet and it might change, and for example there will probably be also a sub-workflows directory and a workflows directory where the pipeline sub-workflows and workflows are also defined.

So that is all I wanted to share with you today. I will be happy to take any questions.

</details>
