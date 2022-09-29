---
title: Software licence
subtitle: Pipelines must open source, released with the MIT licence.
menu:
  main:
    weight: 50
---

All nf-core pipelines must be released with an [MIT licence](https://choosealicense.com/licenses/mit/).

Please do not bundle any other scripts within the workflow, in case they have a different or incompatible licence (for example, in the `bin` directory).
If you need such a script, even a simple one, please release it on bioconda instead and reference it like any other software.
