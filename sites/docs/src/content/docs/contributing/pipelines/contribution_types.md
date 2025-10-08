---
title: Contribution types
subtitle: Guidance for distinguishing between different types of contributors to nf-core pipelines
parentWeight: 10
---

All nf-core pipelines have a `nextflow.config` file, that within which has a 'manifest' listing all contributors to a pipeline.

This section of the manifest have a specific attribute called 'contributions' that is designed to describe what each person has done for the pipeline.

In this page we will describe the nf-core community's guidance on how to use the contribution attribution, to maximise community recognition to all of out pipelines.

## Contribution types

There are three main 'contribution' types that can be specified in this section of the manifest, as described in the [Nextflow documentation](https://nextflow.io/docs/latest/reference/config.html#manifest):

- Author
- Maintainer
- Contributor

We will now explain how nf-core pipeline developers should use these in their pipelines.

### Author

An 'author' is typically someone who was a part of the team that designed and wrote the original version version of the pipeline (i.e., prior and up to version 1), and was therefore the primary maintainer of the pipeline in the past.

The other use of the 'author' designation is for someone who was at one point (but is no longer) one of the lead maintainers of the pipeline.
This means that in the past they ensured the pipeline remained up to date, reviewed PRs, developed a long-term roadmap, added new features, and help coordinate development with other developers and users.

Therefore there can be multiple authors.

### Maintainer

The maintainer designation is someone who is the _currently_ active as a lead developer of the pipeline, and is regularly contributing new code and managing the project.
This means they are the main person or team member of the people who is ensuring the pipeline is up to date, reviewing PRs, is executing a long-term roadmap, adding new features, and currently coordinating with other developers and users.
Maintainers typically have the best current 'overview' of the pipeline.

The maintainer designation does not need to be someone who originally authored the pipeline, but it can be the same.
Once a 'maintainer' steps down as a lead and is no longer regularly or actively contributing or managing the pipeline, they are moved to the 'author' designation (optionally with the 'contributor' designation - see below).

### Contributor

A contributor is someone who only occasionally contributes to the pipeline (whether that is code or docs), and does not have a 'mangement' role.
They do not typically triage and deal with GitHub issues (although they may help contribute bug fixes on an irregular basis).

Contributors can be both current and also in the past.

## Case study

We can use nf-core/mag as an example of distinguishing between the different types of contributions

![](images/contributing/pipelines/mag_contributors_plot.png)

- Hadrien was the original architect and wrote the majority of v1 of the pipeline from 2018, and finished developing v1 in 2020. He is therefore an author.
- Daniel and Sabrina started assisting Hadrien in 2019, and took over as the lead maintainers between 2019-2022 keeping the pipeline up to date and adding new features. During this period they were maintainers.
- Maxime and James started in occasioanlly contributing feautures and bug fixes between 2020-2021, but were not handling the project management. During this period they were just contributors.
- In 2023, Sabrina moved on and stepped down as a project lead. She then moved to an 'author' designation.
- When Sabrina stepped down, James took over with managing the roadmap and triaging bug fixes, he thus changed from 'contributor' to 'maintainer' status.

## Final Notes

It is important to note that a person can be designated one or multiple of the contribution types above!
This should be combined on a case-by-case basis.

Otherwise, as a final point - it is important that no one should ever be _removed_ from this list!
Everyone who has contributed at any point in the history of a pipeline must be acknowledged and recognised, regardless if it was in v1 or v10 of a pipeline.

nf-core is all about community and collaboration.
We hope with this guide we can help strengthen community ownership of all of our pipelines, and give them a long lifetime!
