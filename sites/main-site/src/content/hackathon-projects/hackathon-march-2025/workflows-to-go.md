---
title: Workflows to go - nf-core pipelines download
category: tooling
slack: "https://nfcore.slack.com/archives/C08JHLV9QRJ"
intro_video: "https://youtu.be/5TceBWlltoE"
image: "/assets/images/events/2025/hackathon-march/workflows-to-go.png"
image_alt: A hexagonal photo of a hiking backpack in a rocky landscape, captioned with "nf-core to go" above it.
leaders:
  MatthiasZepper:
    name: Matthias Zepper
    slack: "https://nfcore.slack.com/team/U02ENQANGF9"
---

Let's face it: we are not exactly a light-traveling company! Gathering everything a workflow needs for offline use results in a daunting packing list: the workflow itself, containerized software, configurations, plugins, Nextflow, Java, and more.

Fortunately, `nf-core pipelines download` is here to help with the packing. However, this sub-tool is due for some maintenance, refactoring, and code cleanup. As the final step in a long chain of dependencies, it must constantly adapt to the latest pipeline templates while ideally maintaining backward compatibility with older releases.

Navigating this nearly 2,000 lines codebase isn't easy. But don’t worry — you’ll get a comprehensive introduction from an experienced "mountain" guide to help you scale this mountain of code. And if you’d like, we can even try pair-coding at times.

## Goals

The [full backlog](https://github.com/nf-core/tools/issues?q=is%3Aissue%20state%3Aopen%20download%20label%3Adownload) - peek at your own risk ;-)

1. **Amenable**: Remove hard dependency on Singularity for better MacOS support, rather warn about incomplete container downloads.
2. **Amenable**: Improve container registry error messages: [#2892](https://github.com/nf-core/tools/issues/2892)
3. **Advanced**: Refactor `WorkflowRepo` class and remove `SyncedRepo` superclass: [#2940](https://github.com/nf-core/tools/issues/2940). Get inspiration from [here](https://github.com/aws-samples/amazon-omics-tutorials/blob/main/utils/scripts/nf/__init__.py).
4. **Advanced**: Download plugins with the workflows: [#3344](https://github.com/nf-core/tools/issues/3344)
5. **Ambitious**: Wrap `nextflow inspect` and run on revised `WorkflowRepo`
6. **Astounding**: Move container downloads to a separate tool that is optionally called by `nf-core pipelines download`: [#2408](https://github.com/nf-core/tools/issues/2408)
7. **Astounding**: Support download/set-up of Docker containers as well [see this post](https://github.com/nextflow-io/nextflow/discussions/4708).
8. **Audacious**: Full multi-arch support for Seqera Containers [#3179](https://github.com/nf-core/tools/issues/3179)
9. **Audacious**: Download Galaxy Singularity containers through CVMFS support than via http: [#18972](https://github.com/galaxyproject/galaxy/issues/18972#issuecomment-2404725035)

## Who might like this project?

Feels like cheating, but is totally legal - dodge Nextflow at a Nextflow Hackathon and work with Python instead!

This project is perfect for participants looking to dip their toes into the nf-core ecosystem without diving deep into Nextflow itself! The nf-core tools are implemented in Python, making them an ideal entry point e.g. for scientists who already use the language to analyze their data and want to leverage their skills to contribute to nf-core.

This opportunity is for you if you’re more interested in deploying nf-core pipelines on various infrastructures rather than developing the workflows themselves. While contributing, you'll gain valuable experience on the [DevOps](https://about.gitlab.com/topics/devops/) side of nf-core and learn how to deploy pipelines with confidence on your HPC. Your work will help make life easier for scientists worldwide who need to run pipelines in offline environments!
