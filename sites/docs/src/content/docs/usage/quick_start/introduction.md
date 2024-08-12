---
title: Quick start
subtitle: Get started with Nextflow and nf-core
shortTitle: Introduction
weight: 1
parentWeight: 10
---

## What is nf-core?

nf-core is a community effort to collect a curated set of analysis pipelines built using [Nextflow](https://www.nextflow.io/docs/latest/index.html).

nf-core provides a standardised set of best practices, guidelines, and templates for building and sharing bioinformatics pipelines. Pipelines are designed to be modular, scalable, and portable, allowing researchers to easily adapt and execute them using their own data and compute resources.

One of the key benefits of nf-core is that it promotes open development, testing, and peer review, ensuring that the pipelines are robust, well-documented, and validated against real-world datasets. This helps to increase the reliability and reproducibility of bioinformatics analyses and ultimately enables researchers to accelerate their scientific discoveries.

### Key Features

- **Documentation**
  - nf-core pipelines have extensive documentation covering installation, usage, and description of output files to ensure that you won't be left in the dark.
- **CI Testing**
  - Every time a change is made to the pipeline code, nf-core pipelines use continuous-integration testing to ensure that nothing has broken.
- **Stable Releases**
  - nf-core pipelines use GitHub releases to tag stable versions of the code and software, making pipeline runs totally reproducible.
- **Packaged software**
  - Pipeline dependencies are automatically downloaded and handled using Docker, Singularity, Conda, or other software management tools. There is no need for any software installations.
- **Portable and reproducible**
  - nf-core pipelines follow best practices to ensure maximum portability and reproducibility. The large community makes the pipelines exceptionally well-tested and easy to execute.
- **Cloud-ready**
  - nf-core pipelines are tested on AWS after every major release. You can even browse results live on the website and use outputs for your own benchmarking.

nf-core is published in [Nature Biotechnology: Nat Biotechnol 38, 276–278 (2020)](https://www.nature.com/articles/s41587-020-0439-x). Nature Biotechnology. An updated preprint is available at [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.05.10.592912v1).

## Pipelines

There are currently >110 nf-core pipelines at various stages of development.

Each released pipeline has a dedicated page that includes 6 documentation sections:

- **Introduction**: An introduction and overview of the pipeline
- **Usage**: Descriptions of how to execute the pipeline
- **Parameters**: Grouped pipeline parameters with descriptions
- **Output**: Descriptions and examples of the expected output files
- **Results**: Example output files generated from the full test dataset
- **Releases & Statistics**: pipeline version history and statistics

Each section should be explored to understand what the pipeline does and how it can be configured.

### Join the nf-core community!

The nf-core Slack is the place to get help, share expertise, and connect with other community members.

<p>
    <a class="btn btn-success me-2" href="/join/slack"
        ><i class="fab fa-slack"></i> Join the nf-core Slack</a
    >
</p>

We hope to see you there!

**Important to know**

- **Community conduct**: It’s important that this community is an open, collaborative, and welcoming environment. Read our nf-core [Code of Conduct](https://nf-co.re/code_of_conduct).
- **Channels**: The nf-core Slack has a channel for just about everything! Common topics and each pipeline has its own channel. When you join, you will automatically be signed up to some of our favorite channels. You can leave and join groups at any time.
- **Github invitations**: Post your GitHub username on the [`#github-invitations`](https://nfcore.slack.com/archives/CEB982K2T) channel and we will send you an invitation to join the nf-core GitHub organization.
- **Your Slack profile**: Don’t forget to fill in your profile (especially your GitHub username) so your accounts are linked.
- **Nextflow forum**: There is also a Nextflow Community Forum for asking Nextflow questions. Follow [this link](https://community.seqera.io/) to join.
- **Say hey!**: Introduce yourself in [`#say-hello`](https://nfcore.slack.com/archives/C04F9UNM46M) If you need inspiration, community members often share a few fun notes about their backgrounds, their goals, and some of their favorite things.

There are also several groups that promote diversity and inclusion, e.g., `#nf-women-nb`, please join any that apply to you.

:::note
You’ll probably end up joining lot of channels, so it can be a good idea to organize your sidebar with custom sections.
:::

:::warning
Avoid using `@here` and `@channel` unless absolutely necessary. The `@here` command will notify all currently active channel members. The `@channel` command will send a message to all channel members, regardless of their activity status. If you are responding to a question or comment, try to respond in a thread to help keep conversations organised and together, especially if the channel is busy with other posts.
:::
