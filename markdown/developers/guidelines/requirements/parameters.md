---
title: Standardised parameters
subtitle: Strive to have standardised usage.
menu:
  main:
    weight: 100
---

Where possible pipelines should use the same command line option names as other pipelines for comparable options.
For example, `--input` and `--genome`.

In addition to the names of paramters, they should ideally work in a similar way.
For example, `--input` typically takes a `.csv` sample sheet file (but not always, where not appropriate).

> We're planning to build a tool that lists every paramter used by every pipeline, so that you can check for existing parameters with similar names
> You can track the progress of this feature request here: []()
