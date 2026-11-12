---
title: Static typing in nf-core
category: components
slack: https://nfcore.slack.com/channels/hack-barcelona-2026-static-typing
location: Barcelona
image: "/assets/images/events/2026/hackathon-barcelona/static-typing.jpg"
image_alt: "Meme about a guy complaining that the amount of untyped processes is too damn high."
leaders:
  nvnieuwk:
    name: Nicolas Vannieuwkerke
    slack: https://nfcore.slack.com/team/U03CKGEU3LZ
---

## Goal

Experiment with the new static typing feature in Nextflow processes, write some guidelines for nf-core and prepare the tooling.

## Description

[Static typing in Nextflow processes](https://docs.seqera.io/nextflow/reference/process/inputs-outputs-typed) has been a preview feature for a while now. This will however become a stable feature once Nextflow 26.10.0 releases.
Lots of people in nf-core have already started experimenting using the new typing and the opinions are quite split on how to implement this properly into nf-core. This project will try to compile all ideas and suggestions and create one set of guidelines to rule them all!

Anyone with any opinion or ideas on the implementation of static typing in processes is welcome to join us, no matter your skill level.

We might also handle the implementation of record types if the discussions end sooner than expected, but this is considered to be an optional goal of the project.

## Tasks

- Experiment with static typing in processes and start forming some more opinions on it's implementation in nf-core
- Write guidelines for the community to follow once static typing is officially part of nf-core
- Update nf-core tools to work with statically typed processes
- [OPTIONAL] Repeat the above steps on record types once all other work is done

## Non-goals

We are not planning on updating all modules by the end of the hackathon. This project is aimed towards discussion of the topic and preparing the guidelines and tooling for the offical adoption some time after the hackathon.
