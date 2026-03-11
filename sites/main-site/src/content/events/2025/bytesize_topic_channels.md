---
title: "Bytesize: topic channels in nf-core/modules"
subtitle: Nicolas Vannieuwkerke, Center for Medical Genetics Ghent
type: talk
startDate: "2025-12-02"
startTime: "13:00+01:00"
endDate: "2025-12-02"
endTime: "13:30+01:00"
locations:
  - name: Online
    links:
      - https://kth-se.zoom.us/j/68390542812
---

This week, Nicolas ([@nvnieuwk](https://github.com/nvnieuwk)) is introducing topic channels for nf-core/modules.
A topic channel is a channel that can receive values from many sources implicitly based on a matching topic name.
You can think of it as a channel that is shared across many different processes using the same topic name.
Using a topic channel is as such a convenient way to collect related items from many different sources without explicitly connecting them (e.g. using the mix operator).

Resources:
[Nextflow documentation](https://nextflow.io/docs/latest/reference/channel.html#topic)
