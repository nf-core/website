---
title: Module container conversion party
category: components
slack: https://nfcore.slack.com/archives/CJRH30T6V
location: Barcelona
image: /assets/images/events/2026/hackathon-barcelona/container-conversion.png
image_alt: "Evergreen container ship stuck in the Suez Canal with the label 'nf-core/modules with biocontainers'. Next to it is a small excavator (small compared to the huge ship) digging soil away from the ship to release it. This excavator is labeled as 'hackathon group'."
leaders:
  mribeirodantas:
    name: Marcel Ribeiro-Dantas
    slack: https://nfcore.slack.com/team/U03932BSX1V
  mashehu:
    name: Matthias Hörtenhuber
    slack: https://nfcore.slack.com/team/UQZG2UCF3
---

## Goal

Convert at least one module from BioContainers to Seqera containers

## Description

- New(ish) to nf-core and don't know what to do?
- Trying to get started with something?
- Interested in making nf-core modules more reproducible?
- Want to try out the newest nf-core tools features?

If any of the above is true, then this group is for you.

We will go through [nf-core/modules](https://github.com/nf-core/modules) and swap their containers from BioContainers to Seqera containers using

```bash
nf-core modules container create <module_name>
```

Should be straight forward, and work straight out of the box for most modules, but this unlocks many new possibilities (running pipelines on ARM out-of-the-box, automated tool version bump, world dominance?).

## Tasks

- Select a module from the [nf-core/modules repository](https://github.com/nf-core/modules) that currently uses BioContainers
- [Generate an issue on nf-core/modules](https://github.com/nf-core/modules/issues/new?template=update_module.yml) to claim the module
- Assign yourself to the issue
- Make sure you have the latest `dev` version of nf-core tools installed

```bash
pip install --upgrade --force-reinstall git+https://github.com/nf-core/tools.git@dev
```

- Run `nf-core modules container create <module_name>{:bash}` to generate Wave-built Seqera containers
- Create a pull request to submit the converted module back to nf-core/modules
