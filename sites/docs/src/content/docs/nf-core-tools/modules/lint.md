---
title: Lint a module
subtitle: Check a module against nf-core guidelines
shortTitle: lint
weight: 80
---

Run the `nf-core modules lint` command to check modules in the current working directory (pipeline or nf-core/modules clone) against nf-core guidelines.

- Use the `--all` flag to run linting on all modules found.
- Use `--dir <pipeline_dir>` to specify another directory than the current working directory.
- Use the `--fix` flag to fix all possible linting errors:
  - Write the correct input and output channels and channel structure. This will try to reuse the information existing in your `meta.yml`, and only change the structure of the file if needed. Make sure to check that the added changes are correct and update `description`, `pattern` and `type` values accordingly.
  - Add the missing ontology entries to the input and output entries in `meta.yml`. It will add the proper ontology entries based on the specified `pattern`. Make sure to add any missing ontology URLs.
  - Add a [`bio.tools`](https://bio.tools/) identifier for the tool.

<!-- RICH-CODEX
working_dir: tmp/modules
before_command: sed 's/1.13a/1.10/g' modules/nf-core/multiqc/main.nf > modules/nf-core/multiqc/main.nf.tmp && mv modules/nf-core/multiqc/main.nf.tmp modules/nf-core/multiqc/main.nf
-->

![`nf-core modules lint multiqc`](../../../../assets/images/tools/nf-core-modules-lint.svg)
