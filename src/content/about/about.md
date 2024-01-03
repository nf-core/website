---
title: About nf-core
description: About nf-core - who we are, how we are organised and how we started.
---

# Introduction

The nf-core project is a diverse project spread across many groups.
Here you can read how we organise ourselves, how we are funded and how the nf-core project was started.

Please note that all nf-core community members are expected to adhere to our [code of conduct](/code_of_conduct).

# Financial Support

The vast majority of nf-core development is done as a labour of love, on a voluntary basis.
Most of all, we would like to thank all contributors (and their employers!).
You form the lifeblood of nf-core and we are eternally grateful for your time and efforts.

A number of projects and grants list nf-core as collaborators and also contribute to our community (see [_Projects we are involved with_](/contributors#initiatives)).
If you work with a project that would benefit from an explicit link with nf-core, please let us know.

Finally, we would specifically like to acknowledge and thank the following sponsors who help to support the project:

## Chan Zuckerberg Initiative

<img src="/images/contributors/colour/CZI.svg" alt="Chan Zuckerberg Initiative" class="float-end darkmode-image me-5 mb-5 w-25 ms-3">

The Chan Zuckerberg Initiative (CZI) is a great supporter of scientific open-source software.
We are very grateful to them for supporting nf-core and [Nextflow](https://nextflow.io/) with several grants through their _Essential Open Source Software for Science_ (EOSS) grants:

- [EOSS - second round](https://chanzuckerberg.com/eoss/proposals/nextflow-and-nf-core/)
- [EOSS - fourth round](https://chanzuckerberg.com/eoss/proposals/nextflow-and-nf-core-reproducible-workflows-for-the-scientific-community-cycle-4/)
- [EOSS - Diversity & Inclusion](https://cziscience.medium.com/advancing-diversity-and-inclusion-in-scientific-open-source-eaabe6a5488b)

Amongst other things, the CZI EOSS grant money has enabled the nf-core community to:

- Hire a dedicated coordinator / safety officer (50%)
- Hire a dedicated developer for nf-core (50%)
- Cover operational costs, such as our Slack subscription fees
- Organise events, such as hackathons
- Run organised mentorship programmes
- Hiring local community advocates from geographical areas where we currently lack representation

The CZI grant employs personnel based at the [SciLifeLab National Genomics Infrastructure](https://ngisweden.scilifelab.se/), the [Quantitative Biology Center](http://qbic.life/) in Tübingen Germany, and [Seqera Labs](https://seqera.io/).

<div class="clearfix"></div>

## Seqera Labs

<img src="/images/contributors/colour/seqera.svg" alt="Seqera Labs" class="float-end darkmode-image me-5 mb-5 w-25 ms-3">

Seqera Labs is the leading provider of open source workflow orchestration software needed for data pipeline processing, cloud infrastructure, and secure collaboration.

Seqera is the company behind Nextflow and has supported the nf-core community since its inception.

Specifically, Seqera helps nf-core with:

- Event organisation and costs
- Website hosting
- Employment of several core contributors
- A [Tower Cloud Professional](https://cloud.tower.nf/pricing/) account for launching and managing full-size release tests.
- Providing and maintaining Nextflow! ✨

<div class="clearfix"></div>

## SciLifeLab Data Centre

<img src="/images/contributors/colour/SciLifeLabDC.svg" alt="SciLifeLab Data Centre" class="float-end darkmode-image me-5 mb-5 w-25 ms-3">

The [SciLifeLab Data Centre](https://www.scilifelab.se/data/) supports nf-core with funding from the [SciLifeLab & Wallenberg National Program for Data-Driven Life Science](https://www.scilifelab.se/data-driven/).
This funding covers a full-time position ([@mashehu](https://github.com/mashehu), Matthias Hörtenhuber) to work on maintainance of nf-core framework code. For example, the [nf-core website](https://nf-co.re/) and the [nf-core/tools package](https://github.com/nf-core/tools/).

<div class="clearfix"></div>

## Amazon Web Services

<img src="/images/contributors/colour/aws.svg" alt="Amazon Web Services" class="float-end darkmode-image me-5 mb-5 w-25 ms-3" style="max-width: 200px">

Amazon Web Services (AWS) kindly support nf-core with cloud compute credits to run each nf-core analysis pipeline with full-size benchmark datasets on every release.
You can explore and download these pipeline results under the <em class="mx-2"><i class="fab fa-aws me-2"></i> Results</em> tab on each pipeline page.

AWS also hosts the [AWS-iGenomes](https://registry.opendata.aws/aws-igenomes/) resource on the [Registry of Open Data on AWS](https://registry.opendata.aws/).
This is used by most nf-core pipelines to give free and open access to the reference genomes of over 30 species, by using a simple `--genome` key when running a pipeline.

<div class="clearfix"></div>

## Microsoft Azure

<img src="/images/contributors/colour/azure.svg" alt="Microsoft Azure" class="float-end darkmode-image me-5 mb-5 ms-3">

Microsoft Azure also kindly supports nf-core with cloud compute credits to run each nf-core analysis pipeline with full-size benchmark datasets on every release.
You will soon be able to explore and download these pipeline results on the nf-core website, on each pipeline page.

<div class="clearfix"></div>

# Open Source Support

We also thank the following organisations for supporting nf-core through providing us 'open source plans' of their services:

<div class="row row-cols-1 row-cols-lg-3 my-5">
  <div class="col" style="align-self:center;">
    <a href="https://www.docker.com/" target="_blank">
      <img src="/images/docker-horizontal.png" alt="Docker" class="w-75">
    </a>
  </div>
  <div class="col" style="align-self:center;">
    <a href="https://gitpod.io/" target="_blank">
      <img src="/images/contributors/colour/gitpod.svg" alt="Gitpod" class="darkmode-image w-75">
    </a>
  </div>
  <div class="col" style="align-self:center;">
    <a href="https://hackmd.io/" target="_blank">
      <img src="/images/contributors/colour/hackmd.svg" alt="HackMD" class="darkmode-image w-75">
    </a>
  </div>
</div>

# History of nf-core

The nf-core project came about at the start of 2018. [Phil Ewels](http://phil.ewels.co.uk/) ([@ewels](https://github.com/ewels/)) was the head of the development facility at [NGI Stockholm](https://ngisweden.scilifelab.se/) (National Genomics Infrastructure), part of [SciLifeLab](https://www.scilifelab.se/) in Sweden.

The NGI had been developing analysis pipelines for use with it's genomics data for several years and started using a set of standards for each pipeline created. This helped other people run the pipelines on their own systems; typically Swedish research groups at first, but later on other groups and core genomics facilities too such as [QBIC](http://qbic.life/) in Tübingen.

As the number of users and contributors grew, the pipelines began to outgrow the SciLifeLab and NGI branding. To try to open up the effort into a truly collaborative project, [nf-core](https://github.com/nf-core) was created and all relevant pipelines moved to this new GitHub Organisation.

The early days of nf-core were greatly shaped by [Alex Peltzer](https://apeltzer.github.io/) ([@apeltzer](https://github.com/apeltzer/)), [Sven Fillinger](https://uni-tuebingen.de/en/research/research-infrastructure/quantitative-biology-center-qbic/team0/sven-fillinger/) ([@sven1103](https://github.com/sven1103/)) and [Andreas Wilm](https://andreas-wilm.github.io/) ([@andreas-wilm](https://github.com/andreas-wilm/)).
Without them, the project would not exist.
