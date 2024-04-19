---
title: Subworkflow commands
subtitle: Subworkflows
weight: 0
---

After the launch of nf-core modules, we can provide now also nf-core subworkflows to fully utilize the power of DSL2 modularization.
Subworkflows are chains of multiple module definitions that can be imported into any pipeline.
This allows multiple pipelines to use the same code for a the same tasks, and gives a greater degree of reusability and unit testing.

To allow us to test modules and subworkflows together we put the nf-core DSL2 subworkflows into the `subworkflows` directory of the modules repository is at <https://github.com/nf-core/modules>.
