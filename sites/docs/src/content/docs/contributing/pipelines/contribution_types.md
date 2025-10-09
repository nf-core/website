---
title: Contribution types
subtitle: Guidance for distinguishing between different types of contributors to nf-core pipelines
parentWeight: 10
---

All nf-core pipelines have a `nextflow.config` file, that within which has a 'manifest' listing all contributors to a pipeline.

This section of the manifest have a specific attribute called 'contributions' that is designed to describe what each person has done for the pipeline.

In this page we will describe the nf-core community's guidance on how to use the contribution attribution, to maximise community recognition to all of out pipelines.

## Contribution types

The [manifest]((https://nextflow.io/docs/latest/reference/config.html#manifest)) supports three contribution types:

- `author`
- `maintainer`
- `contributor`

The following sections describe how nf-core attributes contributions using these types.

### Author

An `author` is someone who designed and wrote the original pipeline (prior to and including version 1.0) and served as the primary maintainer at that time.

Someone may also be an `author` if they previously served as a lead `maintainer` of the pipeline. As a lead `maintainer`, they:

- Kept the pipeline up to date
- Reviewed pull requests
- Developed the long-term roadmap
- Added new features
- Coordinated development with other developers and users

A pipeline can have multiple authors.

### Maintainer

A `maintainer` is currently active as a lead developer of the pipeline and regularly contributes code and manages the project. Maintainers:

- Keep the pipeline up to date
- Review pull requests
- Execute the long-term roadmap
- Add new features
- Coordinate with other developers and users

A `maintainer` typically has the best current understanding of the pipeline.
A pipeline can have multiple maintainers.

A `maintainer` doesn't need to be an original author of the pipeline, but can be.

When a `maintainer` steps down and no longer actively contributes to or manages the pipeline, they become an `author`.
You can optionally also add the `contributor` designation (see [Contributor](#contributor)).

### Contributor

A `contributor` occasionally contributes to the pipeline (code, documentation, or other) but doesn't have a management role. Contributors typically do not triage or resolve GitHub issues, but may occasionally contribute features or resolve bugs.

Contributors can be both current and past.

## Case study

The nf-core/mag pipeline demonstrates how to distinguish between the different contribution types.

![](@assets/images//mag_contributors_plot.png)

- Hadrien was the original architect and wrote the majority of v1 of the pipeline from 2018, and finished developing v1 in 2020. He is therefore an author.
- Daniel and Sabrina started assisting Hadrien in 2019, and took over as the lead maintainers between 2019-2022 keeping the pipeline up to date and adding new features. During this period they were maintainers.
- Maxime and James started in occasioanlly contributing feautures and bug fixes between 2020-2021, but were not handling the project management. During this period they were just contributors.
- In 2023, Sabrina moved on and stepped down as a project lead. She then moved to an 'author' designation.
- When Sabrina stepped down, James took over with managing the roadmap and triaging bug fixes, he thus changed from 'contributor' to 'maintainer' status.

## Final notes

A person can have one or more contribution designations.
Evaluate each contributor's role and apply the appropriate types.

:::note{title="Important"}
Never _remove_ anyone from the contributor list. Acknowledge everyone who has contributed to the pipeline at any point in its history, whether in version 1 or version 10.
:::

Community and collaboration are core to nf-core.
Using these contribution types helps strengthen community ownership of pipelines and supports their long-term maintenance!
