---
title: nf-core/tools config builder
category: tooling
slack: "https://nfcore.slack.com/archives/C05UZM8L1C0"
location: Sydney (University of Sydney)
image: "/assets/images/events/2025/hackathon-march/Australian-Biocommons-Logo.png"
image_alt: "Australian BioCommons logo"
leaders:
  mgeaghan:
    name: Michael Geaghan
    slack: "https://nfcore.slack.com/team/U08NCPR1T9N"
---

## Description

Nextflow is a highly configurable workflow language, enabling it to be run on a wide variety of hardware, from local laptops to high performance computing and cloud infrastructure. Configuring Nextflow to optimise it for specific data or hardware requires writing one or more configuration files. To make this process more accessible and user-friendly, we are developing a new nf-core tools utility for quickly generating these config files. Using the Python-based Textual framework, the new config builder utility will present users with a series of questions and fields to fill out in order to build either infrastructure or pipeline module configuration files.

## Goals

- Develop code for reliably writing Nextflow configuration files using user-supplied information.
- Write Pydantic validation code to check user input and prevent invalid configurations.
- Create automations for identifying Nextflow processes and labels available for configuration.
- Testing and improving code for automatic HPC environment detection, including schedulers, container software, and module systems.
- Improving and refactoring code for better performance.

This is an opportunity to contribute to the open-source nf-core tools utility and improve the accessibility and ease-of-use of Nextflow. We welcome contributors of all experience levels.
