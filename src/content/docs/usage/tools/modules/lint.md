---
title: Lint a module
subtitle: Check a module against nf-core guidelines
weight: 80
---

Run the `nf-core modules lint` command to check modules in the current working directory (pipeline or nf-core/modules clone) against nf-core guidelines.

Use the `--all` flag to run linting on all modules found. Use `--dir <pipeline_dir>` to specify another directory than the current working directory.

<!-- RICH-CODEX
working_dir: tmp/modules
before_command: sed 's/1.13a/1.10/g' modules/multiqc/main.nf > modules/multiqc/main.nf.tmp && mv modules/multiqc/main.nf.tmp modules/multiqc/main.nf
-->

![`nf-core modules lint multiqc`](/images/tools/nf-core-modules-lint.svg)
