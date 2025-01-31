---
title: Software license
subtitle: Pipelines must open source, released with the MIT license.
shortTitle: MIT license
menu:
  main:
    weight: 50
---

All nf-core pipelines must be released with an [MIT license](https://choosealicense.com/licenses/mit/).

Please do not bundle any other scripts within the workflow, in case they have a different or incompatible license (for example, in the `bin` directory).
If you need such a script, even a simple one, please release it on bioconda instead and reference it like any other software.

If adding to bioconda is not possible, custom script files written explicitly for the pipeline MUST include a comment at the top with the license information matching that of the pipeline.
See [here](https://github.com/nf-core/mag/blob/61c23fe8914503acabf537bef3682f5b84dc7bc0/bin/combine_tables.py#L3-L4) for an example.
