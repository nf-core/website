---
title: "Bytesize: the `nf-core download` command"
subtitle: Matthias Zepper, National Genomics Infrastructure, Stockholm
type: talk
startDate: "2024-07-09"
startTime: "13:00+02:00"
endDate: "2024-07-09"
endTime: "13:30+02:00"
youtubeEmbed: https://youtu.be/5TceBWlltoE
locations:
  - name: Online
    links:
      - https://youtu.be/5TceBWlltoE
---

Today's bytesize talk can be regarded as a continuation of and supplement to an ongoing series of talks focused on container image management in computational workflows. Over the past weeks, you already heard about the [`nextflow inspect` subcommand to list the container images required by a Nextflow pipeline](https://nf-co.re/events/2024/bytesize_singularity_containers_hpc) and learned how to best use [Singularity](https://nf-co.re/events/2024/bytesize_singularity_containers_hpc) and [Wave](https://nf-co.re/events/2024/bytesize_using_wave) for your containerization needs.

To round off the series, Matthias Zepper ([@MatthiasZepper](https://github.com/MatthiasZepper)) will showcase the use and functionality of the `nf-core download` CLI tool. It has been developed with the hope to streamline running nf-core pipelines in offline environments; a common requirement when processing sensitive genomic information, e.g. for clinical applications. The talk will cover two exemplary scenarios how to set up an offline environment and provide some guidance on the required additional configuration. Time-permitting, you may optionally receive a brief overview of the developmental roadmap and possible future features.
