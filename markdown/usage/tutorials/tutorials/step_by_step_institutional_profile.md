---
title: Step-by-step guide to writing an institutional profile
subtitle: Walkthrough on what you need to set up an nf-core institutional profile
---

## Introduction

nf-core offers a centralised place to store Nextflow configuration profiles that work at an _institutional_ level.

What this means is that you can specify common nextflow pipeline configurations and options that can be shared across all users of that particular institutional cluster.

nf-core offers two level of profile sharing: institutional and pipeline-specific profiles.

- Institutional configuration represent configuration options that apply to all users of _all_ nf-core pipelines. These typically define settings regarding the cluster itself, such as the type of scheduler being used, maximum resource limits and so on.
- Pipeline-specific profiles that represent configuration options that apply to all users of a _specific_ nf-core pipelines. These typically define common parameters all users of the pipeline would use, or customise resource requirements for particular processes of that specific pipeline

This walkthrough will guide you through making an _institutional_-level profile. It will:

- Describe commonly useful information that is worth gathering _before_ writing such a profile
- Go through each component of the profile to describe how this should be written
- Show how to test such a profile before submitting it to nf-core

## Preparation

The first thing you can do is go through the following checklist

### Do you have permission to make the profile public

nf-core institutional profiles are stored publicly on the [nf-core/configs](https://github.com/nf-core/configs/) repository. 

In some cases, system administrators of the cluster at your institution may which to keep certain aspects of the cluster for security reasons. We therefore recommend you check with your sysadmins that you have permission to make such a profile and submit it to nf-core. We recommend you send the link to the repository with one of the configs files to show as an example of the sort of information that would be posted.

### What container engines does your cluster offer

nf-core highly recommends the use of container engines or software environment for running truly reproduble pipelines. This means the actual tools (with correct versions) used within the pipeline are contained in a singular 'image' file.

Therefore, you need to find out what container engines/environments your cluster offers. For nf-core pipelines to work, you need one of the following: `docker`, 

### What are the resource maximums

### Does your cluster user a scheduler

### Do you have to specify 