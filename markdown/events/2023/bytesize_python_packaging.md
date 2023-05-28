---
title: 'Bytesize: Converting Python scripts into packages for PyPI, Bioconda & Biocontainers'
subtitle: Phil Ewels - Seqera Labs
type: talk
start_date: '2023-05-02'
start_time: '13:00 CEST'
end_date: '2023-05-02'
end_time: '13:30 CEST'
youtube_embed: https://www.youtube.com/watch?v=hOuS6mXCwhk
location_url:
  - https://www.youtube.com/watch?v=hOuS6mXCwhk
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: Converting Python scripts into packages for PyPI, Bioconda & Biocontainers

This week, Phil Ewels ([@ewels](https://github.com/ewels/)) will show you how to take a Python script and turn it into a stand-alone command-line tool, ready for distribution via the [Python Package Index](https://pypi.org/) (PyPI).

> You can download a `.zip` file of the "before" and "after" code examples Phil demoed [here](/assets/markdown_assets/events/2023/bytesize-python-packaging/python-packaging.zip).

This is a good thing to do for a few reasons:

- More people can use your scripts - not just within Nextflow
  - This is useful for development, for stand-alone testing
  - It's useful for people using other workflow managers
  - It helps when users are testing a method / debugging with small sample sizes
- It allows scripts to be released under different licenses to the pipeline itself
- Software packaging, that is providing container images with all requirements, is handled automatically

Even if it's a small script that you think no-one will ever use outside of your pipeline, it's easy to do and you don't lose anything 🙂

Once released in PyPI, releases via [Bioconda](https://bioconda.github.io/) are simple (see [Bytesize 40: Software packaging](https://nf-co.re/events/2022/bytesize-40-software-packaging)).
Once in BioConda, software will be available for Conda users, but also Docker + Singularity, via the [BioContainers](https://biocontainers.pro/) project.

<details markdown="1"><summary>Video transcription</summary>
**Note: The content has been edited for reader-friendliness**

[0:01](https://www.youtube.com/watch?v=hOuS6mXCwhk&t=1)

</details>
