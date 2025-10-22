---
title: nf-co2footprint plugin development & integration 🌱
category: tooling
slack: https://nfcore.slack.com/channels/nf-co2footprint
image: /assets/images/events/2025/hackathon-march/globalcore-stripes.png
image_alt: The average annual global temperature over the years 1850-2017, known as the 'warming stripes' figure from the [climate lab book](https://www.climate-lab-book.ac.uk/2018/warming-stripes/) website
leaders:
  skrakau:
    name: Sabrina Krakau
    slack: https://nfcore.slack.com/team/UMLKFJ264
  JosuaCarl:
    name: Josua Carl
    slack: https://nfcore.slack.com/archives/D08JJNCQ4SJ
  nadnein:
    name: Nadja Volkmann
    slack: https://nfcore.slack.com/archives/D08KCG9SRSL
---

## Summary

The nf-co2footprint plugin estimates the CO₂ equivalent emissions of workflow runs. After our first stable release, we are now focusing on integrating the plugin into as many nf-core pipelines as possible to make energy and carbon tracking a standard feature across the community.
In parallel, we are implementing several final refinements and additional features for the next release, and we warmly invite contributions and feedback from the community.
Our goal is to raise individual and collective awareness within the nf-core community about the environmental impact of computational research and to encourage more sustainable software development practices.

## Goals

- Update your configs with parameters for accurate CO₂ emission estimations
- Advance the plugin with respect to its integration into nf-core and cloud platforms

## Tasks

- **⚙️ Add nf-co2footprint to your config files**  
  This can be as easy as adding the plugin ID to your institution’s config, but accuracy may improve by supplying additional parameters.  
  Let’s do it together and discuss how this process could be streamlined for other plugin users.

- **☁️ Improve support for AWS and other cloud providers**  
  This includes several sub-tasks such as:
  - Adding S3 output support
  - Automatically detecting the active AWS region and instance type
  - Extending the CPU TDP reference table with corresponding cloud hardware specifications

- **Align the plugin config syntax with standard Nextflow style**  
  Ensure consistency, readability, and seamless integration with existing configurations.

- **Reduce CO₂ footprint of nf-co2footprint itself**  
  Evaluate and optimize the plugin’s own computational impact to make its operation as efficient and sustainable as possible.

- **nf-core 🤝 nf-co2footprint plugin**  
  As our long-term goal is to make the plugin an integral part of the **nf-core** framework, we are always open to suggestions, feedback, and ideas that help us move in that direction.

## Resources

- [Website](https://nextflow-io.github.io/nf-co2footprint/)
- [Github](https://github.com/nextflow-io/nf-co2footprint)
- [Slack](https://nfcore.slack.com/channels/nf-co2footprint)
