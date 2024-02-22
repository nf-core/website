---
title: "Bytesize: Removing the 'lib/' directory from the template"
subtitle: 'Sayonara lib! Harshil Patel, Seqera'
type: talk
startDate: '2024-02-27'
startTime: '13:00+01:00'
endDate: '2024-02-27'
endTime: '13:30+01:00'
locationURL: https://seqera-io.zoom.us/j/5968645264
---

The latest release of nf-core/tools (v2.13) brings with it some significant changes.
One of the most visible is the removal of the pipeline `lib/` folder.
This directory previously contained groovy code that was imported into workflows to handle common functions
such as handling completion emails, pipeline citations and printing params summaries to the terminal.
Going forward, this code will be included within a shared subworkflow instead, hopefully easing future template syncs.
This should allow for a more modular and flexible codebase, which will allow for more frequent updates and easier maintenance.

Harshil Patel ([@drpatelh](https://github.com/drpatelh)) is giving us an overview of all the changes and how pipeline developers can adapt to the new structure.
